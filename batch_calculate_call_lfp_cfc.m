function cfcResults = batch_calculate_call_lfp_cfc(baseDir,eData,varargin)

pnames = {'chunk_size_s','overlap_size_s','selectFilt','avgTetrodes','sessType','ds_factor','filtBank','fs','avgTrial','n_trial_boot','tRange','crossBat'};
dflts  = {0.5,0.25,[1 2],false,'communication',2,[],2083,true,0,[],false};
[chunk_size_s,overlap_size_s,selectFilt,avgTetrodes,sessionType,ds_factor,filtBank,fs,avgTrial,n_trial_boot,tRange,crossBat] = internal.stats.parseArgs(pnames,dflts,varargin{:});

if isempty(filtBank)
    filtBank = get_HG_theta_filters(fs,ds_factor);
end

switch sessionType
    case 'communication'
        lfp_fnames = dir(fullfile(baseDir,'*LFP_call_trig.mat'));
    case 'operant'
        lfp_fnames = dir(fullfile(baseDir,'*LFP_call_trig_operant.mat'));
end

callTypes = {'calling','listening'};

cfcResults = struct('MIstruct',[],'batNum',[],'sessionType',[],'pair_bat_num',[],...
    'nTrial',[],'call_bat_nums',[],'used_call_IDs',[],'expDates',[],'callType',[]);

cfc_k = 1;
nBat = length(eData.batNums);
[all_session_lfp, all_call_bat_nums, all_callIDs, all_exp_dates, lfp_call_offset] = load_call_lfp_data(eData,lfp_fnames,avgTetrodes,ds_factor,tRange);
t = tic;
for bat_k = 1:nBat
    batNum = eData.batNums{bat_k};
    for call_type_k = 1:length(callTypes)
        if strcmp(callTypes{call_type_k},'calling')
            select_call_idx = strcmp(all_call_bat_nums{bat_k},eData.batNums{bat_k});
        else
            select_call_idx = ~strcmp(all_call_bat_nums{bat_k},eData.batNums{bat_k});
        end
        
        [current_session_lfp, current_call_bat_nums, current_callIDs, current_exp_dates, pair_bat_nums] =...
            get_current_lfp(eData,batNum,crossBat,select_call_idx,all_session_lfp,all_call_bat_nums,all_callIDs,all_exp_dates);
        
        nPair = length(current_session_lfp);
        
        for pair_k = 1:nPair
            [cfcResults(cfc_k).MIstruct, cfcResults(cfc_k).nTrial] = calculate_call_lfp_cfc(current_session_lfp{pair_k},...
                filtBank,lfp_call_offset,'chunk_size_s',chunk_size_s,'overlap_size_s',overlap_size_s,...
                'selectFilt',selectFilt,'avgTrial',avgTrial,'n_trial_boot',n_trial_boot);
            
            cfcResults(cfc_k).call_bat_nums = current_call_bat_nums{pair_k};
            cfcResults(cfc_k).expDates = current_exp_dates{pair_k};
            cfcResults(cfc_k).used_call_IDs = current_callIDs{pair_k};
            cfcResults(cfc_k).pair_bat_num = pair_bat_nums{pair_k};
            
            cfcResults(cfc_k).batNum = batNum;
            cfcResults(cfc_k).callType = callTypes{call_type_k};
            cfcResults(cfc_k).sessionType = sessionType;
            fprintf('%d bats analyzed, %d s elapsed',cfc_k,round(toc(t)));
            cfc_k = cfc_k + 1;
        end
    end
end

end

function [all_call_lfp, all_call_bat_nums, all_callIDs, all_exp_dates, lfp_call_offset] = load_call_lfp_data(eData,lfp_fnames,avgTetrodes,ds_factor,tRange)

lfp_file_strs = arrayfun(@(x) strsplit(x.name,'_'),lfp_fnames,'un',0);
batNums = cellfun(@(x) x{1},lfp_file_strs,'un',0);
expDates = cellfun(@(x) datetime(x{2},'InputFormat','yyyyMMdd'),lfp_file_strs,'un',0);

nBat = length(eData.batNums);
[all_call_lfp, all_call_bat_nums, all_callIDs, all_exp_dates] = deal(cell(1,nBat));
for bat_k = 1:nBat
    batIdx = strcmp(batNums,eData.batNums{bat_k});
    bat_lfp_fnames = lfp_fnames(batIdx);
    bat_exp_dates = expDates(batIdx);
    n_bat_exp = sum(batIdx);
    
    [session_lfp, call_bat_nums, callIDs, current_exp_dates] = deal(cell(1,n_bat_exp));
    activeChannels = eData.activeChannels{bat_k};
    for exp_k = 1:n_bat_exp
        current_lfp_fname = fullfile(bat_lfp_fnames(exp_k).folder,bat_lfp_fnames(exp_k).name);
        [session_lfp{exp_k}, call_bat_nums{exp_k}, callIDs{exp_k}, lfp_call_offset] = prepare_call_lfp_data_for_cfc(current_lfp_fname,activeChannels,'avgTetrode',avgTetrodes,'ds_factor',ds_factor,'tRange',tRange);
        current_exp_dates{exp_k} = repmat(bat_exp_dates(exp_k),1,length(callIDs{exp_k}));
    end
    
    all_call_lfp{bat_k} = cat(3,session_lfp{:});
    all_call_bat_nums{bat_k} = [call_bat_nums{:}];
    all_callIDs{bat_k} = [callIDs{:}];
    all_exp_dates{bat_k} = [current_exp_dates{:}];
end

end

function [current_session_lfp, current_call_bat_nums, current_callIDs, current_exp_dates, pair_bat_nums] = get_current_lfp(eData,batNum,crossBat,select_call_idx,all_session_lfp,all_call_bat_nums,all_callIDs,all_exp_dates)
bat_k = strcmp(eData.batNums,batNum);
if ~crossBat
    current_session_lfp = {all_session_lfp{bat_k}(:,:,select_call_idx)};
    current_call_bat_nums = {all_call_bat_nums{bat_k}(select_call_idx)};
    current_callIDs = {all_callIDs{bat_k}(select_call_idx)};
    current_exp_dates = {all_exp_dates};
    pair_bat_nums = {batNum};
else
    pair_bat_nums = setdiff(eData.batNums,batNum);
    nPair = length(pair_bat_nums);
    used_call_IDs = all_callIDs{bat_k}(select_call_idx);
    [current_session_lfp, current_call_bat_nums, current_callIDs, current_exp_dates] = deal(cell(1,nPair));
    for pair_k = 1:nPair
        pair_bat_idx = strcmp(eData.batNums,pair_bat_nums{pair_k});
        pair_call_idx = cell(1,2);
        [current_callIDs{pair_k},pair_call_idx{:}] = intersect(used_call_IDs,all_callIDs{pair_bat_idx});
        
        combined_session_lfp = {all_session_lfp{bat_k},all_session_lfp{pair_bat_idx}};
        
        for bat_pair_k = 1:2
            combined_session_lfp{bat_pair_k} = combined_session_lfp{bat_pair_k}(:,:,pair_call_idx{bat_pair_k});
        end
        current_session_lfp{pair_k} = combined_session_lfp;
        current_call_bat_nums{pair_k} = all_call_bat_nums{bat_k}(pair_call_idx{1});
        current_exp_dates{pair_k} = all_exp_dates{bat_k}(pair_call_idx{1});
    end
end

end

function filtBank = get_HG_theta_filters(fs,ds_factor)

fs_ds = fs/ds_factor;
stop_band_freqs = [1 8; 90 210];
band_pass_freqs = [3 6; 100 200];
attenuationDb = 30*ones(size(band_pass_freqs,1),2);
pass_ripple = ones(1,size(band_pass_freqs,1));
filtType = 'bandpassfir';
filtBank = generate_filter_bank(filtType,band_pass_freqs,stop_band_freqs,attenuationDb,pass_ripple,fs_ds);

end

