function  [MIstruct,nTrial,t] = calculate_call_lfp_cfc(session_lfp,batNum,call_bat_nums,filtBank,include_call_flag,lfp_call_offset,varargin)

pnames = {'chunk_size_s','selectFilt','n_phase_bins','nSurr','randtype','fs','ds_factor'};
dflts  = {[],[],18,0,2,2083,2};
[chunk_size_s,selectFilt,nbins,nSurr,randtype,fs,ds_factor] = internal.stats.parseArgs(pnames,dflts,varargin{:});

minCalls = 10;

fs_ds = round(fs/ds_factor);

nT = size(session_lfp,2);
t = linspace(-lfp_call_offset,lfp_call_offset,nT);
if isempty(chunk_size_s)
    t_idx = 1:nT;
else
    chunk_size_samples = round(fs_ds*chunk_size_s);
    t_idx = slidingWin(nT,chunk_size_samples,0);
end
t = mean(t(t_idx),2)';

nChunks = size(t_idx,1);
MIstruct = cell(1,nChunks);

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
if nTrial < minCalls
    [MIstruct{:}] = deal(struct('MI',[],'MIp',[],'MIsurrmi',[],'MIsurr',[],'mean_amps',[],'ninds',[]));
    return
end
for t_k = 1:nChunks
    current_lfp = session_lfp(:,t_idx(t_k,:),callIdx);
    MIstruct{t_k} = calculate_cfc(current_lfp,filtBank,'selectFilt',selectFilt,...
        'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
end


end