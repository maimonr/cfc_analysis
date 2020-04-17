function  MIstruct = calculate_cfc(current_lfp,filtBank,varargin)

pnames = {'selectFilt','n_phase_bins','nSurr','randtype','filt_chunk_size','parFlag'};
dflts  = {[],18,0,2,200,false};
[selectFilt,nbins,nSurr,randtype,filt_chunk_size,parFlag] = internal.stats.parseArgs(pnames,dflts,varargin{:});

miParams = struct('nbins',nbins,'nSurr',nSurr,'randtype',randtype);

if iscell(current_lfp)
    multi_bat_flag = true;
    assert(length(current_lfp) == 2)
else
    multi_bat_flag = false;
end

nFilt = length(filtBank);
if isempty(selectFilt)
    if nFilt > 50
        disp('This will result in very many filter combinations, are you sure you want to proceed?')
        keyboard
    end
    
    if ~multi_bat_flag
        selectFilt = combnk(1:nFilt,2);
    else
        [A,B] = meshgrid(1:nFilt,nFilt + (1:nFilt));
        c = cat(2,A',B');
        selectFilt = reshape(c,[],2);
    end
end

nSelectFilt = size(selectFilt,1);

filteredHilb = get_filtered_hilbert_lfp(current_lfp,filtBank);

if nSelectFilt > filt_chunk_size
    MIstruct = chunk_MI_calculation_by_filter(filteredHilb,selectFilt,filt_chunk_size,miParams,parFlag);
else
    hilbAmp = abs(filteredHilb(:,:,:,selectFilt(:,2)));
    hilbPhase = angle(filteredHilb(:,:,:,selectFilt(:,1)));
    MIstruct = wrap_get_mi(hilbAmp,hilbPhase,miParams,parFlag);
end
    
MIstruct = reshape(vertcat(MIstruct{:}),size(MIstruct));

end

function filteredHilb = get_filtered_hilbert_lfp(current_lfp,filtBank)

if ~iscell(current_lfp)
    current_lfp = {current_lfp};
    nBat = 1;
else
    nBat = length(current_lfp);
end

nFilt = length(filtBank);
filteredHilb = cell(1,nBat);

for bat_k = 1:nBat
    
    nChannel = size(current_lfp{bat_k},1);
    nSamp = size(current_lfp{bat_k},2);
    nTrial = size(current_lfp{bat_k},3);
    
    filteredHilb{bat_k} = nan([nChannel nSamp nTrial nFilt]);
    for ch_k = 1:nChannel
        current_channel_lfp = squeeze(current_lfp{bat_k}(ch_k,:,:));
        for filt_k = 1:nFilt
            if ~any(isnan(current_channel_lfp),'all')
                filteredPower = filtfilt(filtBank{filt_k},current_channel_lfp);
                filteredHilb{bat_k}(ch_k,:,:,filt_k) = hilbert(squeeze(filteredPower));
            end
        end
    end
end

filteredHilb = cat(4,filteredHilb{:});

end

function MIstruct = wrap_get_mi(hilbAmp,hilbPhase,miParams,parFlag)

nChannel = size(hilbAmp,1);
chunkSize = size(hilbAmp,4);
MIstruct = cell(nChannel,chunkSize);

if parFlag
    parfor filt_k = 1:chunkSize
        current_hilb_amp = hilbAmp(:,:,:,filt_k);
        current_hilb_phase = hilbPhase(:,:,:,filt_k);
        for ch_k = 1:nChannel
            miInput_amp = reshape(squeeze(current_hilb_amp(ch_k,:,:)),1,[])';
            miInput_phase = reshape(squeeze(current_hilb_phase(ch_k,:,:)),1,[])';
            MIstruct{ch_k,filt_k} = get_mi(miInput_phase,miInput_amp,miParams.nbins,miParams.nSurr,miParams.randtype); %#ok<PFBNS>
        end
    end
else
    for filt_k = 1:chunkSize
        current_hilb_amp = hilbAmp(:,:,:,filt_k);
        current_hilb_phase = hilbPhase(:,:,:,filt_k);
        for ch_k = 1:nChannel
            miInput_amp = reshape(squeeze(current_hilb_amp(ch_k,:,:)),1,[])';
            miInput_phase = reshape(squeeze(current_hilb_phase(ch_k,:,:)),1,[])';
            MIstruct{ch_k,filt_k} = get_mi(miInput_phase,miInput_amp,miParams.nbins,miParams.nSurr,miParams.randtype);
        end
    end
end

end

function MIstruct = chunk_MI_calculation_by_filter(filteredHilb,selectFilt,filt_chunk_size,miParams,parFlag)

nChannel = size(filteredHilb,1);
nSelectFilt = size(selectFilt,1);

filtChunks = slidingWin(nSelectFilt,filt_chunk_size,0);
MIstruct = cell(nChannel,nSelectFilt);

n_filt_chunk = size(filtChunks,1);
cum_filt_k = 1;

for f_chunk_k = 1:n_filt_chunk
    
    filt_chunk_idx = selectFilt(filtChunks(f_chunk_k,:),:);
    
    if f_chunk_k == n_filt_chunk
        remaining_chunk_idx = selectFilt(filtChunks(end)+1:nSelectFilt,:);
        filt_chunk_idx = [filt_chunk_idx; remaining_chunk_idx]; %#ok<AGROW>
    end
    
    hilb_amp_chunk = abs(filteredHilb(:,:,:,filt_chunk_idx(:,2)));
    hilb_phase_chunk = angle(filteredHilb(:,:,:,filt_chunk_idx(:,1)));
    
    chunkSize = size(filt_chunk_idx,1);
    
    MIstruct(:,cum_filt_k:(cum_filt_k+chunkSize-1)) = wrap_get_mi(hilb_amp_chunk,hilb_phase_chunk,miParams,parFlag);
    
    cum_filt_k = cum_filt_k + chunkSize;
end


end