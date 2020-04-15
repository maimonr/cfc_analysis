function [MICorr, MICorr_boot, bat_pair_nums] = get_session_MI_corr(all_session_cfc_results,expDate)

nBoot = 1e3;

[ts,MI,batNums] = get_aligned_session_MI(all_session_cfc_results,expDate);

L = length(ts);

nBat = length(MI);
batPairs = combnk(1:nBat,2);
nPair = size(batPairs,1);

bat_pair_nums = batNums(batPairs);

[MICorr, MICorr_boot] = deal(cell(1,nPair));
parfor pair_k = 1:nPair
    X1 = MI{batPairs(pair_k,1)};
    X2 = MI{batPairs(pair_k,2)};
    r = corr(X1,X2);
    MICorr{pair_k} = r(:);
    MICorr_boot{pair_k} = nan(nBoot,length(MICorr{pair_k}));
    for boot_k = 1:nBoot
        shift = randi(L,1);
        r = corr(circshift(X1,shift),X2);
        MICorr_boot{pair_k}(boot_k,:) = r(:);
    end
end

end