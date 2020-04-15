function  [MIstruct,nTrial,tRange] = calculate_call_lfp_cfc(session_lfp,call_bat_nums,exp_type,filtBank,include_call_flag,varargin)

pnames = {'tRange','selectFilt','concatChannels','n_phase_bins','nSurr','randtype'};
dflts  = {[],[],false,18,0,2};
[tRange,selectFilt,concatChannels,nbins,nSurr,randtype] = internal.stats.parseArgs(pnames,dflts,varargin{:});

minCalls = 10;
eData = ephysData(exp_type);
nBat = length(eData.batNums);
ds_factor = 2;
fs_orig = 2083;
fs_ds = round(fs_orig/ds_factor);
lfp_all_avg_ds = cell(1,nBat);
for bat_k = 1:nBat
    if isempty(session_lfp{bat_k}) || ndims(session_lfp{bat_k}) < 3 %#ok<ISMAT>
        continue
    end
    sz = size(session_lfp{bat_k});
    lfp_all_avg_ds{bat_k} = nan(sz(1),round(sz(2)/ds_factor),sz(3));
    for ch_k = 1:sz(1)
        lfp_tmp = squeeze(session_lfp{bat_k}(ch_k,:,:));
        lfp_all_avg_ds{bat_k}(ch_k,:,:) = resample(lfp_tmp,1,ds_factor);
    end
end

%%

if isempty(tRange)
    nT = 1;
else
    nT = size(tRange,2);
end

nTrial = zeros(1,nBat);
MIstruct = cell(nT,nBat);

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
        [MIstruct{:,bat_k}] = deal(struct('MI',[],'MIp',[],'MIsurrmi',[],'MIsurr',[],'mean_amps',[],'ninds',[]));
        continue
    end
    for t_k = 1:nT
        if isempty(tRange)
            t_idx = true(1,size(lfp_all_avg_ds{bat_k},2));   
        else
            t_idx = round(tRange(1,t_k)*fs_ds):round(tRange(2,t_k)*fs_ds);
        end
        
        current_lfp = lfp_all_avg_ds{bat_k}(:,t_idx,idx);
        MIstruct{t_k,bat_k} = calculate_cfc(current_lfp,filtBank,'selectFilt',selectFilt,...
            'concatChannels',concatChannels,'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
    end
end

end