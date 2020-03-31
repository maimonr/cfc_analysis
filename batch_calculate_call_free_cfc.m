function [call_free_cfc_results,success,err] = batch_calculate_call_free_cfc(baseDir,eData,cData,filtBank,outDir)

total_dirs = 1;
t = tic;
success = [];
exp_dirs = dir(fullfile(baseDir,'*20*'));

results_fname = fullfile(outDir,[datestr(datetime,'yyyymmdd') '_call_free_cfc_results.mat']);

all_exp_dirs = cell(1,length(exp_dirs));
for k = 1:length(exp_dirs)
    all_exp_dirs{k} = dir(fullfile(exp_dirs(k).folder,exp_dirs(k).name,'lfpformat','*LFP.mat'));
end

fs = 2083;
notch_filter_60Hz=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5,'DesignMethod','butter','SampleRate',fs);
notch_filter_120Hz=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',119.5,'HalfPowerFrequency2',120.5,'DesignMethod','butter','SampleRate',fs);
notchFilters = {notch_filter_60Hz,notch_filter_120Hz};

call_free_cfc_results = struct('MIstruct',[],'batNum',[],'expDate',[]);

lastProgress = 0;
for k = 1:length(exp_dirs)
    exp_dir = fullfile(exp_dirs(k).folder,exp_dirs(k).name);
    for bat_k = 1:length(eData.batNums)
        lfp_data_fname = dir(fullfile(exp_dir,'lfpformat',[eData.batNums{bat_k} '*LFP.mat']));
        if ~isempty(lfp_data_fname)
            expDate = strsplit(lfp_data_fname.name,'_');
            expDate = datetime(expDate{2},'InputFormat','yyyyMMdd');
            lfp_data_fname = fullfile(lfp_data_fname.folder,lfp_data_fname.name);
            usedChannels = eData.activeChannels{bat_k};
            try
                lfpData = load(lfp_data_fname);
                call_free_cfc_results(total_dirs).MIstruct = calculate_call_free_cfc(lfpData,cData,expDate,filtBank,notchFilters,usedChannels);
                call_free_cfc_results(total_dirs).expDate = expDate;
                call_free_cfc_results(total_dirs).batNum = eData.batNums{bat_k};
            catch err
                errs{total_dirs} = err;
                success(total_dirs) = false;
            end
            success(total_dirs) = true;
            
            progress = 100*(total_dirs/length(all_exp_dirs));
            elapsed_time = round(toc(t));
            
            if mod(progress,10) < mod(lastProgress,10)
                save(results_fname,'call_free_cfc_results')
                fprintf('%d %% of current bat''s directories  processed\n',round(progress));
                fprintf('%d total directories processed, %d s elapsed\n',total_dirs,elapsed_time);
            end
            total_dirs = total_dirs + 1;
            lastProgress = progress;
        end
    end
end