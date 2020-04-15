function MIstruct = calculate_all_session_cfc(lfpData,filtBank,notchFilters,varargin)

pnames = {'avgTetrodes','selectFilt','n_phase_bins','nSurr','randtype','ds_factor','call_offset_length','cfcType','chunk_size_s','fs','cData','expDate','artifact_nStd_factor'};
dflts  = {false,[1 2],18,0,2,2,5,'callFree',2,2083,[],[],5};
[avgTetrodes,selectFilt,nbins,nSurr,randtype,ds_factor,call_offset_length,cfcType,chunk_size_s,fs,cData,expDate,artifact_nStd_factor] = internal.stats.parseArgs(pnames,dflts,varargin{:});

if any(strcmp(cfcType,{'callFree','callOnly'}))
    callPos = cData('expDay',expDate).callPos;
    callPos = callPos(:,1);
    sessionIdx = lfpData.timestamps >= min(0,min(callPos));
else
    sessionIdx = lfpData.timestamps >= 0;
end
timestamps = lfpData.timestamps(sessionIdx);
usedChannels = lfpData.active_channels;

lfpData = lfpData.lfpData(:,sessionIdx);

nChannel = size(lfpData,1);
nT = size(lfpData,2);
lfp_ds = nan(nChannel,round(nT/ds_factor));

for ch_k = 1:nChannel
    try
        lfp_ds(ch_k,:) = resample(lfpData(ch_k,:),1,ds_factor);
    catch err
        if strcmp(err.identifier,'MATLAB:subsassigndimmismatch') && ch_k == 1
            current_lfp_ds = resample(lfpData(ch_k,:),1,ds_factor);
            lfp_ds = nan(nChannel,length(current_lfp_ds));
            lfp_ds(ch_k,:) = current_lfp_ds;
        end
    end
end

lfpData = lfp_ds;
clear lfp_ds
nT = size(lfpData,2);
timestamps = downsample(timestamps,ds_factor);
fs_ds = fs/ds_factor;

if any(strcmp(cfcType,{'callFree','callOnly'}))
    
    call_free_idx = true(1,length(timestamps));
    for k = 1:length(callPos)
        [~,callIdx] = inRange(timestamps,[callPos(k)-call_offset_length callPos(k)+call_offset_length]);
        call_free_idx = call_free_idx & ~callIdx;
    end
    
end

if avgTetrodes
    nTt = 4;
    nChannel_per_tt = 4;
    tetrodeChannels = reshape(0:(nTt*nChannel_per_tt)-1,nTt,[]);
    ttLFP = nan(nTt,size(lfpData,2));
    for tt_k = 1:nTt
        used_channel_idx = ismember(usedChannels,tetrodeChannels(:,tt_k));
        current_lfp = squeeze(nanmean(lfpData(used_channel_idx,:),1));
        ttLFP(tt_k,:) = current_lfp;
    end
    lfpData = ttLFP;
    nChannel = nTt;
end

for ch_k = 1:nChannel
    for notch_filt_k = 1:length(notchFilters)
        lfpData(ch_k,:) = filtfilt(notchFilters{notch_filt_k},lfpData(ch_k,:));
    end
end

switch cfcType
    
    case 'wholeSession'
        
         MIstruct = calculate_cfc(lfpData,filtBank,'selectFilt',selectFilt,...
            'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
    
    case 'callFree'
        
        lfpData = lfpData(:,call_free_idx);
        
        MIstruct = calculate_cfc(lfpData,filtBank,'selectFilt',selectFilt,...
            'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
        
    case 'callOnly'
        
        lfpData = lfpData(:,~call_free_idx);
        
        MIstruct = calculate_cfc(lfpData,filtBank,'selectFilt',selectFilt,...
            'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
        
    case 'slidingWin'
        
        chunkSize = round(chunk_size_s * fs_ds);
        timeChunks = slidingWin(nT,chunkSize,0);
        nChunks = size(timeChunks,1);
        n_chunk_samp = size(timeChunks,2);
        nFilt = size(selectFilt,1);
        
        chunk_timestamps = mean(timestamps(timeChunks),2);
        
        mu = mean(lfpData,2);
        sigma = std(lfpData,[],2);
        
        MI = nan(nChunks,nChannel,nFilt);
        n_artifact_times = zeros(nChunks,nChannel);
    
        current_csc = zeros(nChannel,n_chunk_samp,nChunks);
        for chunk_k = 1:nChunks
            current_csc(:,:,chunk_k) = lfpData(:,timeChunks(chunk_k,:));
        end
        
        for chunk_k = 1:nChunks
            chunk_lfp_data = current_csc(:,:,chunk_k);
            n_artifact_times(chunk_k,:) = sum(abs(chunk_lfp_data - mu) > artifact_nStd_factor*sigma,2);
            
            current_MIstruct = calculate_cfc(chunk_lfp_data,filtBank,'selectFilt',selectFilt,...
                'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
            MI(chunk_k,:,:) = arrayfun(@(x) x.MI,current_MIstruct);
        end
        MIstruct = struct('MI',MI,'n_artifact_times',n_artifact_times,'timestamps',chunk_timestamps);
end

end