function cfcResults = batch_calculate_cross_bat_call_cfc(baseDir,eData,varargin)

pnames = {'chunk_size_s','overlap_size_s','selectFilt','avgTetrodes','sessType','included_call_type','ds_factor','filtBank','fs','avgTrial','nBoot'};
dflts  = {0.5,0.25,[1 2],true,'communication','all',2,[],2083,true,0};
[chunk_size_s,overlap_size_s,selectFilt,avgTetrodes,sessionType,included_call_type,ds_factor,filtBank,fs,avgTrial,nBoot] = internal.stats.parseArgs(pnames,dflts,varargin{:});

if isempty(filtBank)
    filtBank = get_HG_theta_filters(fs,ds_factor);
end

switch sessionType
    case 'communication'
        lfp_fnames = dir(fullfile(baseDir,'*LFP_call_trig.mat'));
    case 'operant'
        lfp_fnames = dir(fullfile(baseDir,'*LFP_call_trig_operant.mat'));
end

n_lfp_files = length(lfp_fnames);
lfp_file_strs = arrayfun(@(x) strsplit(x.name,'_'),lfp_fnames,'un',0);
batNums = cellfun(@(x) x{1},lfp_file_strs,'un',0);
expDates = cellfun(@(x) datetime(x{2},'InputFormat','yyyyMMdd'),lfp_file_strs);

cfcResults = cell(1,n_lfp_files);

for exp_k = 1:n_lfp_files
    
    current_lfp_fname = fullfile(lfp_fnames(exp_k).folder,lfp_fnames(exp_k).name);
    activeChannels = eData.activeChannels{strcmp(batNums{exp_k},eData.batNums)};
    [session_lfp, call_bat_nums, callIDs, lfp_call_offset] = prepare_call_lfp_data_for_cfc(current_lfp_fname,avgTetrodes,activeChannels,ds_factor);
    
    expIdx = expDates == expDates(exp_k) & ~strcmp(batNums,batNums{exp_k});
    nPairs = sum(expIdx);
    pair_bat_nums = batNums(expIdx)';
    
    cfcResults_tmp = struct('MIstruct',[],'expDate',expDates(exp_k),...
        'batNum',batNums(exp_k),'pair_bat_num',pair_bat_nums,...
        'sessionType',sessionType,'included_call_type',included_call_type,...
        'nTrial',[],'call_bat_nums',[],'used_call_IDs',[]);

    pair_lfp_fnames = lfp_fnames(expIdx);
    pair_call_idx = cell(1,2);
    
    for pair_lfp_k = 1:nPairs
        current_lfp_fname = fullfile(pair_lfp_fnames(pair_lfp_k).folder,pair_lfp_fnames(pair_lfp_k).name);
        activeChannels = eData.activeChannels{strcmp(pair_bat_nums{pair_lfp_k},eData.batNums)};
        [pair_session_lfp, ~, pair_callIDs] = prepare_call_lfp_data_for_cfc(current_lfp_fname,avgTetrodes,activeChannels,ds_factor);
        
        [~,pair_call_idx{:}] = intersect(callIDs,pair_callIDs);
        combined_session_lfp = {session_lfp,pair_session_lfp};
        
        for pair_k = 1:2
            combined_session_lfp{pair_k} = combined_session_lfp{pair_k}(:,:,pair_call_idx{pair_k});
        end
        pair_lfp_sizes = cellfun(@size,combined_session_lfp,'un',0);
        assert(isequal(pair_lfp_sizes{:}))
        
        current_call_bat_nums = call_bat_nums(pair_call_idx{1});
        current_callIDs = callIDs(pair_call_idx{1});
        
        [cfcResults_tmp(pair_lfp_k).MIstruct, cfcResults_tmp(pair_lfp_k).nTrial, ~, used_call_idx] = calculate_call_lfp_cfc(combined_session_lfp,...
            batNums{exp_k},current_call_bat_nums,filtBank,included_call_type,lfp_call_offset,...
            'chunk_size_s',chunk_size_s,'overlap_size_s',overlap_size_s,'selectFilt',selectFilt,...
            'avgTrial',avgTrial);
        
        cfcResults_tmp(pair_lfp_k).call_bat_nums = current_call_bat_nums(used_call_idx);
        cfcResults_tmp(pair_lfp_k).used_call_IDs = current_callIDs(used_call_idx);
        if nBoot > 1
            MIboot = cell(1,nBoot);
            for boot_k = 1:nBoot
                trial_shuffled_lfp = get_trial_shuffled_lfp(combined_session_lfp{2},current_call_bat_nums,batNums{exp_k});
                combined_session_lfp_shuffle = {session_lfp,trial_shuffled_lfp};
                MI_struct_shuffle = calculate_call_lfp_cfc(combined_session_lfp_shuffle,...
                    batNums{exp_k},current_call_bat_nums,filtBank,included_call_type,lfp_call_offset,...
                    'chunk_size_s',chunk_size_s,'overlap_size_s',overlap_size_s,'selectFilt',selectFilt,...
                    'avgTrial',avgTrial);
                MIboot{boot_k} = cellfun(@(MI) [MI.MI],MI_struct_shuffle,'un',0);
            end
            cfcResults_tmp(pair_lfp_k).MIboot = MIboot;
        end
    end
    
    cfcResults{exp_k} = cfcResults_tmp;
    
end

end

function filtBank = get_HG_theta_filters(fs,ds_factor)

fs_ds = fs/ds_factor;
attenuationDb = [15 15; 20 20];
stop_band_freqs = [1 12; 70 160];
band_pass_freqs = [4 8; 80 150];
filtType = 'bandpassfir';
filtBank = generate_filter_bank(filtType,band_pass_freqs,stop_band_freqs,attenuationDb,[2 2],fs_ds);

end

function trial_shuffled_lfp = get_trial_shuffled_lfp(session_lfp,call_bat_nums,batNum)

callingIdx = find(strcmp(call_bat_nums,batNum));
listeningIdx = find(~strcmp(call_bat_nums,batNum));

trial_shuffled_lfp = session_lfp;

trial_shuffled_lfp(:,:,callingIdx) = trial_shuffled_lfp(:,:,callingIdx(randperm(length(callingIdx))));
trial_shuffled_lfp(:,:,listeningIdx) = trial_shuffled_lfp(:,:,listeningIdx(randperm(length(listeningIdx))));
end


