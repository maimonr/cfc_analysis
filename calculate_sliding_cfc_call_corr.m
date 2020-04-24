function [cfc_pairwise_R, sliding_win_t] = calculate_sliding_cfc_call_corr(cfcResults,expType,corr_win_s,corr_overlap_s)

all_exp_dates = unique([cfcResults.expDate]);
nExp = length(all_exp_dates);
batNums = unique({cfcResults.batNum});
nBat = length(batNums);
bat_pair_idx = combvec(1:nBat,1:nBat)';
bat_pair_idx = bat_pair_idx(~(bat_pair_idx(:,1) == bat_pair_idx(:,2)),:);
nPair = size(bat_pair_idx,1);

ds_factor = 2;
fs_ds = round(2083/ds_factor);
chunk_size_s = 0.5;
overlap_size_s = 0.25;

minCalls = 10;

switch expType
    case 'adult_operant'
        lfp_call_offset = 4;
        nT = 8333;
    case 'adult'
        lfp_call_offset = 3;
        nT = 6250;
end
    
t = linspace(-lfp_call_offset,lfp_call_offset,nT);
chunk_size_samples = round(fs_ds*chunk_size_s);
overlap_size_samples = round(fs_ds*overlap_size_s);
sliding_win_idx = slidingWin(nT,chunk_size_samples,overlap_size_samples);
t = mean(t(sliding_win_idx),2);

sliding_corr_win_size = round(corr_win_s/overlap_size_s);
sliding_corr_win_overlap = round(corr_overlap_s/overlap_size_s);
sliding_corr_idx = slidingWin(length(t),sliding_corr_win_size,sliding_corr_win_overlap);
nWin = size(sliding_corr_idx,1);

sliding_win_t = mean(t(sliding_corr_idx),2);

cfc_pairwise_R = nan(nWin,nExp,nPair);

parfor win_k = 1:nWin
    call_t_idx = sliding_corr_idx(win_k,:);
    for exp_k = 1:nExp
        current_exp_date = all_exp_dates(exp_k);
        for pair_k = 1:nPair
            pair_idx = bat_pair_idx(pair_k,:);
            [bat_exp_idx,current_MI_struct,current_used_call_IDs] = deal(cell(1,2));
            for bat_k = 1:length(pair_idx)
                bat_exp_idx{bat_k} = [cfcResults.expDate] == current_exp_date & strcmp({cfcResults.batNum},batNums{pair_idx(bat_k)});
                if ~any(bat_exp_idx{bat_k})
                    continue
                end
                current_MI_struct{bat_k} = cfcResults(bat_exp_idx{bat_k}).MIstruct;
                current_used_call_IDs{bat_k} = cfcResults(bat_exp_idx{bat_k}).used_call_IDs;
            end
            if any(cellfun(@isempty,current_MI_struct))
                continue
            end
            
            current_MI_mat = cellfun(@(batMI) cellfun(@(MIstruct) median([MIstruct.MI]),batMI),current_MI_struct,'un',0);
            used_call_idx = cell(1,2);
            [~,used_call_idx{:}] = intersect(current_used_call_IDs{:});
            
            used_call_bat_nums = cfcResults(bat_exp_idx{2}).call_bat_nums(used_call_idx{2});
            
            calling_bat_num = batNums(pair_idx(2));
            callIdx = strcmp(used_call_bat_nums,calling_bat_num);
            current_MI_mat = cellfun(@(MI,idx) MI(idx,:),current_MI_mat,used_call_idx,'un',0);
            if sum(callIdx) < minCalls
                continue
            end
            
            callMI = cellfun(@(MI) reshape(MI(callIdx,call_t_idx),1,[]),current_MI_mat,'un',0);
            
            R = corrcoef(callMI{:});
            cfc_pairwise_R(win_k,exp_k,pair_k) = R(2);
            
        end
    end
end