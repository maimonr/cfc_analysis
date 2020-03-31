function MIstruct = calculate_call_free_cfc(lfpData,cData,expDate,filtBank,notchFilters,usedChannels)

nTt = 4;
nChannel_per_tt = 4;

nbins = 24;
nSurr = 0;
randtype = 2;

callPos = cData('expDay',expDate).callPos;
callPos = callPos(:,1);
sessionIdx = lfpData.timestamps >= min(0,min(callPos));
timestamps = lfpData.timestamps(sessionIdx);
lfpData = lfpData.lfpData(:,sessionIdx);
call_free_idx = true(1,length(timestamps));
call_offset_length = 5;
for k = 1:length(callPos)
    [~,callIdx] = inRange(timestamps,[callPos(k)-call_offset_length callPos(k)+call_offset_length]);
    call_free_idx = call_free_idx & ~callIdx;
end
call_free_lfp = lfpData(:,call_free_idx);
clear lfpData
tetrodeChannels = reshape(0:(nTt*nChannel_per_tt)-1,nTt,[]);
session_lfp = zeros(nTt,sum(call_free_idx));

for ch_k = 1:size(call_free_lfp,1)
    current_lfp = call_free_lfp(ch_k,:);
    for notch_filt_k = 1:length(notchFilters)
        current_lfp = filtfilt(notchFilters{notch_filt_k},current_lfp);
    end
    call_free_lfp(ch_k,:) = current_lfp;
end

for tt_k = 1:nTt
    used_channel_idx = ismember(usedChannels,tetrodeChannels(:,tt_k));
    current_lfp = squeeze(nanmean(call_free_lfp(used_channel_idx,:),1));
    session_lfp(tt_k,:) = current_lfp;
end
clear call_free_lfp

sz = size(session_lfp);
lfp_all_avg_ds = nan(sz(1),round(sz(2)/2));
ds_factor = 2;

for tt_k = 1:nTt
    lfp_tmp = squeeze(session_lfp(tt_k,:));
    lfp_all_avg_ds(tt_k,:) = resample(lfp_tmp,1,ds_factor);
end
clear session_lfp

nFilt = length(filtBank);
filteredHilb = nan([size(lfp_all_avg_ds) nFilt]);
for filt_k = 1:nFilt
    for tt_k = 1:nTt
        current_channel_lfp = squeeze(lfp_all_avg_ds(tt_k,:));
        if ~any(isnan(current_channel_lfp),'all')
            filteredPower = filtfilt(filtBank{filt_k},current_channel_lfp);
            filteredHilb(tt_k,:,filt_k) = hilbert(squeeze(filteredPower));
        end
    end
end

[miInput_amp,miInput_phase] = deal(nan(size(filteredHilb,2),nFilt,nTt));
for filt_k = 1:nFilt
    for tt_k = 1:nTt
        miInput_amp(:,filt_k,tt_k) = abs(squeeze(filteredHilb(tt_k,:,filt_k)));
        miInput_phase(:,filt_k,tt_k) = angle(squeeze(filteredHilb(tt_k,:,filt_k)));
    end
end

MIstruct = get_mi(miInput_phase,miInput_amp,nbins,nSurr,randtype);

end