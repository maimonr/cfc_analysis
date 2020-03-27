function [session_lfp, call_bat_nums] = prepare_data_for_cfc(baseDir,exp_type,session_type,expDate)

eData = ephysData(exp_type);
nBat = length(eData.batNums);

if strcmp(session_type,'communication')
    fName_str = '_LFP_call_trig.mat';
elseif strcmp(session_type,'operant')
    fName_str = '_LFP_call_trig_operant.mat';
end
call_trig_csc = struct('call_trig_csc',[],'batNums',[]);
for bat_k = 1:nBat
    lfpFname = fullfile(baseDir,'lfp_data',[eData.batNums{bat_k} '_' datestr(expDate,'yyyymmdd') fName_str]);
    if exist(lfpFname,'file')
        call_trig_csc(bat_k) = load(lfpFname,'call_trig_csc','batNums');
    else
        call_trig_csc(bat_k) = struct('call_trig_csc',[],'batNums',[]);
    end
end

nTt = 4;
nChannel_per_tt = 4;
session_lfp = cell(1,nBat);

if all(cellfun(@isempty,{call_trig_csc.call_trig_csc}))
    call_bat_nums = cell(1,nBat);
    return
end

for bat_k = 1:nBat
    nT = size(call_trig_csc(bat_k).call_trig_csc,1);
    nTrial = size(call_trig_csc(bat_k).call_trig_csc,2);
    if nTrial < 1
        continue
    end
    session_lfp{bat_k} = zeros(nTt,nT,nTrial);
    usedChannels = eData.activeChannels{bat_k};
    tetrodeChannels = reshape(0:(nTt*nChannel_per_tt)-1,nTt,[]);
    
    for tt_k = 1:nTt
        used_channel_idx = ismember(usedChannels,tetrodeChannels(:,tt_k));
        current_lfp = squeeze(nanmean(call_trig_csc(bat_k).call_trig_csc(:,:,used_channel_idx),3));
        session_lfp{bat_k}(tt_k,:,:) = current_lfp;
    end
end
call_bat_nums = {call_trig_csc.batNums};
end