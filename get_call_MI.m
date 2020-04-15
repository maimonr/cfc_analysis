function [callMI, non_call_MI, nCall] = get_call_MI(baseDir,all_session_cfc_results,expDate)

[ts,MI] = get_aligned_session_MI(all_session_cfc_results,expDate);

[call_t, callIdx, nonCallIdx] = get_session_call_t(baseDir,expDate,ts);

if any(isnan([call_t, callIdx, nonCallIdx]))
    [callMI, non_call_MI, nCall] = deal(NaN);
    return
end

nBat = length(MI);
[callMI, non_call_MI] = deal(cell(1,nBat));
for bat_k = 1:nBat
    callMI{bat_k} = nanmedian(MI{bat_k}(callIdx,:),1);
    non_call_MI{bat_k} = nanmedian(MI{bat_k}(nonCallIdx,:),1);
end
nCall = length(call_t);
end