function [ts,MI,batNums] = get_aligned_session_MI(all_session_cfc_results,expDate)

exp_ks = [all_session_cfc_results.expDate] == expDate;
session_cfc_results = all_session_cfc_results(exp_ks);
batNums = {session_cfc_results.batNum};

ts_orig = arrayfun(@(x) x.MIstruct.timestamps,session_cfc_results,'un',0);
tsRange = [max(cellfun(@min,ts_orig)) min(cellfun(@max,ts_orig))];
ts = cellfun(@(x) inRange(x,tsRange),ts_orig,'un',0);
L = min(cellfun(@length,ts));
ts = cellfun(@(x) x(1:L),ts,'un',0);

assert(length(unique(cellfun(@length,ts))) == 1)
assert(range(cellfun(@range,ts)) < min(cellfun(@(ts) mean(diff(ts)),ts)))

used_ts_idx = cellfun(@(ts,ts_orig) ismember(ts_orig,ts),ts,ts_orig,'un',0);

MI = arrayfun(@(x) x.MIstruct.MI,session_cfc_results,'un',0);
MI = cellfun(@(mi,idx) mi(idx,:),MI,used_ts_idx,'un',0);

ts = mean([ts{:}],2);

end