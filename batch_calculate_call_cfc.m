function cfcResults = batch_calculate_call_cfc(expType)
f = waitbar(0,'Initializing');
band_pass_freqs = [4 8; 80 150];
stop_band_freqs = [2 10; 70 160];
attenuationDb = 30*ones(size(band_pass_freqs,1),2);
pass_ripple = ones(1,size(band_pass_freqs,1));
fs_orig = 2083;
ds_factor = 2;
fs_ds = round(fs_orig/ds_factor);
filtType = 'bandpassfir';
if ~exist('filtBank','var')
    filtBank = generate_filter_bank(filtType,band_pass_freqs,stop_band_freqs,attenuationDb,pass_ripple,fs_ds);
end
%%
avgTetrode = false;
selectFilt = [1 2];
concatChannels = true;
switch expType
    case 'adult_operant'
        baseDir = 'E:\ephys\adult_operant_recording\';
        sessionTypes = {'communication','operant'};
%         tRange = [0.25:0.25:5.5; 2.25:0.25:7.5];
        tRange = [1.5; 6.5];
    case 'adult'
        baseDir = 'E:\ephys\adult_recording\';
        sessionTypes = {'communication'};
%         tRange = [0.25:0.25:3.75; 2.25:0.25:5.75];
        tRange = [0.5; 5.5];
end
timerT = tic;
T = readtable(fullfile(baseDir,'documents','recording_logs.csv'));
expDates = unique(T{logical(T.usable),1});
all_exp_k = 1;
cfcResults = struct('MIStruct',[],'expDate',[],'sessionType',[],'included_call_type',[],'nTrial',[]);
for session_k = 1:length(sessionTypes)
    for exp_k = 1:length(expDates)
        [session_lfp, call_bat_nums] = prepare_call_lfp_data_for_cfc(baseDir,expType,sessionTypes{session_k},expDates(exp_k),avgTetrode);
        
        if strcmp(sessionTypes,'operant')
            include_call_flags = {'all'};
        else
            include_call_flags = {'calling','listening'};
        end
        
        for include_k = 1:length(include_call_flags)
            [cfcResults(all_exp_k).MIstruct, cfcResults(all_exp_k).nTrial] = calculate_call_lfp_cfc(session_lfp,call_bat_nums,...
                expType,filtBank,include_call_flags{include_k},'tRange',tRange,'selectFilt',selectFilt,'concatChannels',concatChannels);
            cfcResults(all_exp_k).expDate = expDates(exp_k);
            cfcResults(all_exp_k).sessionType = sessionTypes{session_k};
            cfcResults(all_exp_k).included_call_type = include_call_flags{include_k};
            currentT = round(toc(timerT));
            waitbar(all_exp_k/(length(expDates)*2*length(sessionTypes)),f,['Calculating CFC ' num2str(currentT) 's elapsed'])
            all_exp_k = all_exp_k + 1;
        end
    end
end
close(f)

end

