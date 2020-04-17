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
lfpData = lfpData.lfpData(:,sessionIdx);
usedChannels = lfpData.active_channels;

[lfpData, timestamps] = get_ds_lfp_data(lfpData,timestamps,ds_factor);

if any(strcmp(cfcType,{'callFree','callOnly'}))
    call_free_idx = get_call_free_timestamp_idx(timestamps,callPos,call_offset_length);
end

if avgTetrodes
    lfpData = average_lfp_by_tt(lfpData,usedChannels);
end

nChannel = size(lfpData,1);

for ch_k = 1:nChannel
    for notch_filt_k = 1:length(notchFilters)
        lfpData(ch_k,:) = filtfilt(notchFilters{notch_filt_k},lfpData(ch_k,:));
    end
end

if strcmp(cfcType,'slidingWin')
    
    MIstruct = calculate_sliding_win_MI(lfpData,fs,ds_factor,chunk_size_s,selectFilt,timestamps,artifact_nStd_factor,filtBank,nbins,nSurr,randtype);
    
else
    switch cfcType
        
        case 'wholeSession'
            % do nothing
        case 'callFree'
            lfpData = lfpData(:,call_free_idx);
        case 'callOnly'
            lfpData = lfpData(:,~call_free_idx);
    end
    
    MIstruct = calculate_cfc(lfpData,filtBank,'selectFilt',selectFilt,...
        'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
end

end

function call_free_idx = get_call_free_timestamp_idx(timestamps,callPos,call_offset_length)

call_free_idx = true(1,length(timestamps));
for k = 1:length(callPos)
    [~,callIdx] = inRange(timestamps,[callPos(k)-call_offset_length callPos(k)+call_offset_length]);
    call_free_idx = call_free_idx & ~callIdx;
end

end

function [lfp_ds, timestamps_ds] = get_ds_lfp_data(lfpData,timestamps,ds_factor)

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

timestamps_ds = downsample(timestamps,ds_factor);
end

function lfpData = average_lfp_by_tt(lfpData,usedChannels)

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

end

function MIstruct = calculate_sliding_win_MI(lfpData,fs,ds_factor,chunk_size_s,selectFilt,timestamps,artifact_nStd_factor,filtBank,nbins,nSurr,randtype)

nChannel = size(lfpData,1);
nT = size(lfpData,2);

fs_ds = fs/ds_factor;
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