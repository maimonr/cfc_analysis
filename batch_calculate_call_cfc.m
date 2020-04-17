function cfcResults = batch_calculate_call_cfc(baseDir,eData,varargin)

pnames = {'chunk_size_s','selectFilt','avgTetrodes','sessType','included_call_type','ds_factor','filtBank','fs'};
dflts  = {2,[1 2],false,'communication','all',2,[],2083};
[chunk_size_s,selectFilt,avgTetrodes,sessionType,included_call_type,ds_factor,filtBank,fs] = internal.stats.parseArgs(pnames,dflts,varargin{:});

f = waitbar(0,'Initializing');

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
expDates = cellfun(@(x) datetime(x{2},'InputFormat','yyyyMMdd'),lfp_file_strs,'un',0);

cfcResults = struct('MIstruct',[],'expDate',expDates,'batNum',batNums,'sessionType',sessionType,'included_call_type',included_call_type,'nTrial',[]);

for exp_k = 1:n_lfp_files
    
    current_lfp_fname = fullfile(lfp_fnames(exp_k).folder,lfp_fnames(exp_k).name);
    activeChannels = eData.activeChannels{strcmp(batNums{exp_k},eData.batNums)};
    [session_lfp, call_bat_nums, lfp_call_offset] = prepare_call_lfp_data_for_cfc(current_lfp_fname,avgTetrodes,activeChannels,ds_factor);
    
    [cfcResults(exp_k).MIstruct, cfcResults(exp_k).nTrial] = calculate_call_lfp_cfc(session_lfp,batNums{exp_k},call_bat_nums,...
        filtBank,included_call_type,lfp_call_offset,'chunk_size_s',chunk_size_s,'selectFilt',selectFilt);
    waitbar(exp_k/n_lfp_files,f,'Calculating CFC');
end

close(f)

end

function filtBank = get_HG_theta_filters(fs,ds_factor)

fs_ds = fs/ds_factor;
attenuationDb = [15 15; 20 20];
stop_band_freqs = [1 12; 70 160];
band_pass_freqs = [4 8; 80 150];
filtType = 'bandpassfir';
filtBank = generate_filter_bank(filtType,band_pass_freqs,stop_band_freqs,attenuationDb,[2 2],fs_ds);

end
