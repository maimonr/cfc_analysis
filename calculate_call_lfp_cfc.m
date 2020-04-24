function  [MIstruct,nTrial,t,callIdx] = calculate_call_lfp_cfc(session_lfp,batNum,call_bat_nums,filtBank,include_call_flag,lfp_call_offset,varargin)

pnames = {'chunk_size_s','overlap_size_s','selectFilt','n_phase_bins','nSurr','randtype','fs','ds_factor','avgTrial'};
dflts  = {[],0,[],18,0,2,2083,2,true};
[chunk_size_s,overlap_size_s,selectFilt,nbins,nSurr,randtype,fs,ds_factor,avgTrial] = internal.stats.parseArgs(pnames,dflts,varargin{:});

minCalls = 10;

fs_ds = round(fs/ds_factor);

if iscell(session_lfp)
    nT = size(session_lfp{1},2);
else
    nT = size(session_lfp,2);
    session_lfp = {session_lfp};
end
nBat = length(session_lfp);

t = linspace(-lfp_call_offset,lfp_call_offset,nT);
if isempty(chunk_size_s)
    t_idx = 1:nT;
else
    chunk_size_samples = round(fs_ds*chunk_size_s);
    overlap_size_samples = round(fs_ds*overlap_size_s);
    t_idx = slidingWin(nT,chunk_size_samples,overlap_size_samples);
end
t = mean(t(t_idx),2)';

nChunks = size(t_idx,1);

switch include_call_flag
    case 'calling'
        callIdx = strcmp(call_bat_nums,batNum);
    case 'listening'
        if ~isempty(call_bat_nums)
            callIdx = ~strcmp(call_bat_nums,batNum);
        else
            callIdx = [];
        end
    case 'all'
        callIdx = true(1,length(call_bat_nums));
end

nTrial = sum(callIdx);

if avgTrial
    used_trial_idx = {callIdx};
else
    used_trial_idx = num2cell(find(callIdx));
end


n_used_trial = length(used_trial_idx);
MIstruct = cell(n_used_trial,nChunks);

if nTrial < minCalls
    [MIstruct{:}] = deal(struct('MI',[],'MIp',[],'MIsurrmi',[],'MIsurr',[],'mean_amps',[],'ninds',[]));
    return
end

for trial_k = 1:n_used_trial
    for t_k = 1:nChunks
        current_lfp = cell(1,nBat);
        for bat_k = 1:nBat
            current_lfp{bat_k} = session_lfp{bat_k}(:,t_idx(t_k,:),used_trial_idx{trial_k});
        end
        MIstruct{trial_k,t_k} = calculate_cfc(current_lfp,filtBank,'selectFilt',selectFilt,...
            'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
    end
end


end