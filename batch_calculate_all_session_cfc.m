function [all_session_cfc_results,success,errs] = batch_calculate_all_session_cfc(baseDir,filtBank,outDir,varargin)

pnames = {'cfcType','cData','chunk_size_s','selectFilt','avgTetrodes','ds_factor','expDate','fs'};
dflts  = {'slidingWin',[],2,[1 2],false,2,[],2083};
[cfcType,cData,chunk_size_s,selectFilt,avgTetrodes,ds_factor,expDate,fs] = internal.stats.parseArgs(pnames,dflts,varargin{:});

cfc_calculation_inputs = cell(1,2*length(dflts));
cfc_calculation_inputs(1:2:end-1) = pnames;
cfc_calculation_inputs(2:2:end) = deal({cfcType,cData,chunk_size_s,selectFilt,avgTetrodes,ds_factor,expDate});

t = tic;

if contains(baseDir,'lfp_data')
    lfp_fnames = dir(fullfile(baseDir,'*LFP.mat'));
else
    expDirs = dir(fullfile(baseDir,'*20*'));
    lfp_fnames = cell(1,length(expDirs));
    
    for exp_k = 1:length(expDirs)
        lfp_fnames{exp_k} = dir(fullfile(expDirs(exp_k).folder,expDirs(exp_k).name,'lfpformat','*LFP.mat'));
    end
    lfp_fnames = vertcat(lfp_fnames{:});
end

[all_session_cfc_results, lfp_fnames, results_fname, n_lfp_files, start_file_idx] = preallocate_cfc_results(lfp_fnames,outDir);

notch_filter_60Hz=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5,'DesignMethod','butter','SampleRate',fs);
notch_filter_120Hz=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',119.5,'HalfPowerFrequency2',120.5,'DesignMethod','butter','SampleRate',fs);
notchFilters = {notch_filter_60Hz,notch_filter_120Hz};

errs = cell(1,n_lfp_files);
success = false(1,n_lfp_files);

lastProgress = 0;
for file_k = start_file_idx + (1:n_lfp_files)
    
    lfp_data_fname = fullfile(lfp_fnames(file_k).folder,lfp_fnames(file_k).name);
    try
        lfpData = load(lfp_data_fname);        
        all_session_cfc_results(file_k).MIstruct = calculate_all_session_cfc(lfpData,filtBank,notchFilters,cfc_calculation_inputs{:});
        success(file_k) = true;
    catch err
        errs{file_k} = err;
    end
    lastProgress = update_progress_and_save(results_fname,lfp_fnames,all_session_cfc_results,file_k,lastProgress,t);
end

end

function [all_session_cfc_results, lfp_fnames, results_fname, n_lfp_files, start_file_idx] = preallocate_cfc_results(lfp_fnames,outDir)

n_lfp_files = length(lfp_fnames);
lfp_file_strs = arrayfun(@(x) strsplit(x.name,'_'),lfp_fnames,'un',0);
batNums = cellfun(@(x) x{1},lfp_file_strs,'un',0);
expDates = cellfun(@(x) datetime(x{2},'InputFormat','yyyyMMdd'),lfp_file_strs,'un',0);

if isfile(outDir)
    results_fname = outDir;
    all_session_cfc_results = load(results_fname);
    all_session_cfc_results = all_session_cfc_results.all_session_cfc_results;
    completedIdx = arrayfun(@(x) ~isempty(x.MIstruct),all_session_cfc_results);
    cfcProgress = all_session_cfc_results(completedIdx);
    cfcProgress = struct2table(rmfield(cfcProgress,'MIstruct'));
    
    to_be_processed_idx = ~ismember(table(batNums,[expDates{:}]','VariableNames',{'batNum','expDate'}),cfcProgress);
    lfp_fnames = lfp_fnames(to_be_processed_idx);
    batNums = batNums(to_be_processed_idx);
    expDates = expDates(to_be_processed_idx);
    all_session_cfc_results = vertcat(all_session_cfc_results(completedIdx), struct('MIstruct',[],'batNum',batNums,'expDate',expDates));
    
    n_lfp_files = sum(to_be_processed_idx);
    start_file_idx = sum(completedIdx);
else
    all_session_cfc_results = struct('MIstruct',[],'batNum',batNums,'expDate',expDates);
    results_fname = fullfile(outDir,[datestr(datetime,'yyyymmdd') '_call_free_cfc_results.mat']);
    start_file_idx = 0;
end

end

function lastProgress = update_progress_and_save(results_fname,lfp_fnames,all_session_cfc_results,file_k,lastProgress,t)

progress = 100*(file_k/length(lfp_fnames));
elapsed_time = round(toc(t));

if mod(progress,10) < mod(lastProgress,10)
    save(results_fname,'all_session_cfc_results')
    fprintf('%d %% of current bat''s directories  processed\n',round(progress));
    fprintf('%d total directories processed, %d s elapsed\n',file_k,elapsed_time);
end
lastProgress = progress;

end