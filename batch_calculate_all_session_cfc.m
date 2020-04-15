function [all_session_cfc_results,success,errs] = batch_calculate_all_session_cfc(baseDir,filtBank,outDir,varargin)

pnames = {'cfcType','cData','chunk_size_s','selectFilt','avgTetrodes','ds_factor','expDate'};
dflts  = {'slidingWin',[],2,[1 2],false,2,[]};
[cfcType,cData,chunk_size_s,selectFilt,avgTetrodes,ds_factor,expDate] = internal.stats.parseArgs(pnames,dflts,varargin{:});

cfc_calculation_inputs = cell(1,2*length(dflts));
cfc_calculation_inputs(1:2:end-1) = pnames;
cfc_calculation_inputs(2:2:end) = deal({cfcType,cData,chunk_size_s,selectFilt,avgTetrodes,ds_factor,expDate});

t = tic;

results_fname = fullfile(outDir,[datestr(datetime,'yyyymmdd') '_call_free_cfc_results.mat']);

lfp_fnames = dir(fullfile(baseDir,'*LFP.mat'));
n_lfp_files = length(lfp_fnames);
lfp_file_strs = arrayfun(@(x) strsplit(x.name,'_'),lfp_fnames,'un',0);

batNums = cellfun(@(x) x{1},lfp_file_strs,'un',0);
expDates = cellfun(@(x) datetime(x{2},'InputFormat','yyyyMMdd'),lfp_file_strs,'un',0);

fs = 2083;
notch_filter_60Hz=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5,'DesignMethod','butter','SampleRate',fs);
notch_filter_120Hz=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',119.5,'HalfPowerFrequency2',120.5,'DesignMethod','butter','SampleRate',fs);
notchFilters = {notch_filter_60Hz,notch_filter_120Hz};

all_session_cfc_results = struct('MIstruct',[],'batNum',batNums,'expDate',expDates);
errs = cell(1,n_lfp_files);
success = false(1,n_lfp_files);

lastProgress = 0;
for file_k = 1:n_lfp_files
    
    lfp_data_fname = fullfile(lfp_fnames(file_k).folder,lfp_fnames(file_k).name);
    try
        lfpData = load(lfp_data_fname);        
        all_session_cfc_results(file_k).MIstruct = calculate_all_session_cfc(lfpData,filtBank,notchFilters,cfc_calculation_inputs{:});
    catch err
        errs{file_k} = err;
        success(file_k) = false;
    end
    success(file_k) = true;
    
    progress = 100*(file_k/length(lfp_fnames));
    elapsed_time = round(toc(t));
    
    if mod(progress,10) < mod(lastProgress,10)
        save(results_fname,'all_session_cfc_results')
        fprintf('%d %% of current bat''s directories  processed\n',round(progress));
        fprintf('%d total directories processed, %d s elapsed\n',file_k,elapsed_time);
    end
    lastProgress = progress;
end

end