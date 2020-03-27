function  [MIstruct,nTrial,tRange] = calculate_cfc(session_lfp,call_bat_nums,exp_type,filtBank,include_call_flag,tRange)

nTt = 4;
minCalls = 10;
eData = ephysData(exp_type);
nBat = length(eData.batNums);
ds_factor = 2;
fs_orig = 2083;
fs_ds = round(fs_orig/ds_factor);
lfp_all_avg_ds = cell(1,nBat);
for bat_k = 1:nBat
    if isempty(session_lfp{bat_k}) || ndims(session_lfp{bat_k}) < 3
        continue
    end
    sz = size(session_lfp{bat_k});
    lfp_all_avg_ds{bat_k} = nan(sz(1),round(sz(2)/ds_factor),sz(3));
    for tt_k = 1:nTt
        lfp_tmp = squeeze(session_lfp{bat_k}(tt_k,:,:));
        lfp_all_avg_ds{bat_k}(tt_k,:,:) = resample(lfp_tmp,1,ds_factor);
    end
end

%%

nbins = 24;
nSurr = 100;
randtype = 2;
nT = size(tRange,2);
nFilt = length(filtBank);
nTrial = zeros(1,nBat);
MIstruct = cell(nT,nBat);
% MIstruct = struct('MI',[],'MIp',[],'MIsurrmi',[],'MIsurr',[],'mean_amps',[],'ninds',[]);
for t_k = 1:nT
    for bat_k = 1:nBat
        switch include_call_flag
            case 'calling'
                idx = strcmp(call_bat_nums{bat_k},eData.batNums{bat_k});
            case 'listening'
                if ~isempty(call_bat_nums{bat_k})
                    idx = ~strcmp(call_bat_nums{bat_k},eData.batNums{bat_k});
                else
                    idx = [];
                end
            case 'all'
                idx = true(1,length(call_bat_nums{bat_k}));
        end
        
        nTrial(bat_k) = sum(idx);
        if nTrial(bat_k) < minCalls
            MIstruct{t_k,bat_k} = struct('MI',[],'MIp',[],'MIsurrmi',[],'MIsurr',[],'mean_amps',[],'ninds',[]);
            continue
        end
        t_idx = round(tRange(1,t_k)*fs_ds):round(tRange(2,t_k)*fs_ds);
        X = lfp_all_avg_ds{bat_k}(:,t_idx,idx);
        nSamp = size(X,2);
        
        filteredHilb = nan([size(X) nFilt]);
        for filt_k = 1:nFilt
            for tt_k = 1:nTt
                current_channel_lfp = squeeze(X(tt_k,:,:));
                if ~any(isnan(current_channel_lfp),'all')
                    filteredPower = filtfilt(filtBank{filt_k},current_channel_lfp);
                    filteredHilb(tt_k,:,:,filt_k) = hilbert(squeeze(filteredPower));
                end
            end
        end
        
        [miInput_amp,miInput_phase] = deal(nan(nTrial(bat_k)*nSamp,nFilt,nTt));
        for filt_k = 1:nFilt
            for tt_k = 1:nTt
                miInput_amp(:,filt_k,tt_k) = reshape(abs(squeeze(filteredHilb(tt_k,:,:,filt_k))),1,[]);
                miInput_phase(:,filt_k,tt_k) = reshape(angle(squeeze(filteredHilb(tt_k,:,:,filt_k))),1,[]);
            end
        end
        MIstruct{t_k,bat_k} = get_mi(miInput_phase,miInput_amp,nbins,nSurr,randtype);
    end
end

end