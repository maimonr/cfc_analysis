function  MIstruct = calculate_cfc(current_lfp,filtBank,varargin)

pnames = {'selectFilt','n_phase_bins','nSurr','randtype'};
dflts  = {[],18,0,2};
[selectFilt,nbins,nSurr,randtype] = internal.stats.parseArgs(pnames,dflts,varargin{:});

nFilt = length(filtBank);
if isempty(selectFilt)
    selectFilt = repmat([1:nFilt]',1,2);
    nSelectFilt = nFilt;
else
    nSelectFilt = size(selectFilt,1);
end

nChannel = size(current_lfp,1);
nSamp = size(current_lfp,2);
nTrial = size(current_lfp,3);

filteredHilb = nan([nChannel nSamp nTrial nFilt]);
for ch_k = 1:nChannel
    current_channel_lfp = squeeze(current_lfp(ch_k,:,:));
    for filt_k = 1:nFilt
        if ~any(isnan(current_channel_lfp),'all')
            filteredPower = filtfilt(filtBank{filt_k},current_channel_lfp);
            filteredHilb(ch_k,:,:,filt_k) = hilbert(squeeze(filteredPower));
        end
    end
end


MIstruct = cell(nChannel,nSelectFilt);

%%

filt_chunk_size = 200;
filtChunks = slidingWin(nSelectFilt,filt_chunk_size,0);
n_filt_chunk = size(filtChunks,1);
tic;
cum_filt_k = 0;
for f_chunk_k = 1:n_filt_chunk
    
    filt_chunk_idx = selectFilt(filtChunks(f_chunk_k,:),:);
    
    if f_chunk_k == n_filt_chunk
        remaining_chunk_idx = selectFilt(filtChunks(end)+1:nSelectFilt,:);
        filt_chunk_idx = [filt_chunk_idx; remaining_chunk_idx]; %#ok<AGROW>
    end
    
    hilb_amp_chunk = abs(filteredHilb(:,:,:,filt_chunk_idx(:,2)));
    hilb_phase_chunk = angle(filteredHilb(:,:,:,filt_chunk_idx(:,1)));
    
    chunkSize = size(filt_chunk_idx,1);
    parfor filt_k = 1:chunkSize
        current_hilb_amp = hilb_amp_chunk(:,:,:,filt_k);
        current_hilb_phase = hilb_phase_chunk(:,:,:,filt_k);
        for ch_k = 1:nChannel
            miInput_amp = reshape(squeeze(current_hilb_amp(ch_k,:,:)),1,[])';
            miInput_phase = reshape(squeeze(current_hilb_phase(ch_k,:,:)),1,[])';
            MIstruct{ch_k,cum_filt_k + filt_k} = get_mi(miInput_phase,miInput_amp,nbins,nSurr,randtype);
        end
    end
    cum_filt_k = cum_filt_k + chunkSize;
end

t = toc;

MIstruct = reshape(vertcat(MIstruct{:}),size(MIstruct));
%%
end