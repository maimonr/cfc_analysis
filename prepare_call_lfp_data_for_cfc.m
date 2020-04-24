function [session_lfp_ds, call_bat_nums, used_call_IDs, lfp_call_offset] = prepare_call_lfp_data_for_cfc(current_lfp_fname,avgTetrode,activeChannels,ds_factor)

call_trig_csc = load(current_lfp_fname,'call_trig_csc','batNums','lfp_call_offset','used_call_IDs');
call_bat_nums = call_trig_csc.batNums;
used_call_IDs = call_trig_csc.used_call_IDs;
lfp_call_offset = call_trig_csc.lfp_call_offset;

call_trig_csc = call_trig_csc.call_trig_csc;

if avgTetrode
    nChannel = 4;
    nChannel_per_tt = 4;
else
    nChannel = size(call_trig_csc,3);
end

nT = size(call_trig_csc,1);
nTrial = size(call_trig_csc,2);

if nTrial < 1
    [session_lfp_ds, call_bat_nums] = deal(NaN);
    return
end

session_lfp = zeros(nChannel,nT,nTrial);

if avgTetrode
    tetrodeChannels = reshape(0:(nChannel*nChannel_per_tt)-1,nChannel,[]);
    for tt_k = 1:nChannel
        used_channel_idx = ismember(activeChannels,tetrodeChannels(:,tt_k));
        current_lfp = squeeze(nanmean(call_trig_csc(:,:,used_channel_idx),3));
        session_lfp(tt_k,:,:) = current_lfp;
    end
else
    session_lfp = permute(call_trig_csc,[3 1 2]);
    if nChannel ~= length(activeChannels) && nChannel == 16
        session_lfp = session_lfp(activeChannels,:,:);
    end
end

session_lfp_ds = nan(nChannel,round(nT/ds_factor),nTrial);
for ch_k = 1:nChannel
    lfp_tmp = squeeze(session_lfp(ch_k,:,:));
    session_lfp_ds(ch_k,:,:) = resample(lfp_tmp,1,ds_factor);
end

end