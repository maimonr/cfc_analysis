function  [MIstruct,nTrial,t] = calculate_call_lfp_cfc(session_lfp,filtBank,lfp_call_offset,varargin)

pnames = {'chunk_size_s','overlap_size_s','selectFilt','n_phase_bins','nSurr','randtype','fs','ds_factor','avgTrial','n_trial_boot','minCalls'};
dflts  = {[],0,[],18,0,2,2083,2,true,0,10};
[chunk_size_s,overlap_size_s,selectFilt,nbins,nSurr,randtype,fs,ds_factor,avgTrial,n_trial_boot,minCalls] = internal.stats.parseArgs(pnames,dflts,varargin{:});

assert((n_trial_boot > 0 && avgTrial) || n_trial_boot == 0)

fs_ds = round(fs/ds_factor);

if iscell(session_lfp)
    nT = size(session_lfp{1},2);
    nTrial = size(session_lfp{1},3);
else
    nT = size(session_lfp,2);
    nTrial = size(session_lfp,3);
    session_lfp = {session_lfp};
end

t = linspace(-lfp_call_offset,lfp_call_offset,nT);
if isempty(chunk_size_s)
    t_idx = 1:nT;
else
    chunk_size_samples = round(fs_ds*chunk_size_s);
    overlap_size_samples = round(fs_ds*overlap_size_s);
    t_idx = slidingWin(nT,chunk_size_samples,overlap_size_samples);
end
t = mean(t(t_idx),2)';
nChunks = size(t_idx,1);

callIdx = true(1,nTrial);

if avgTrial
    n_used_trial = 1;
    used_trial_idx = {callIdx};
else
    n_used_trial = nTrial;
    used_trial_idx = num2cell(find(callIdx));
end

MIstruct = cell(n_used_trial,nChunks);

if nTrial < minCalls
    [MIstruct{:}] = deal(struct('MI',[],'MIp',[],'MIsurrmi',[],'MIsurr',[],'mean_amps',[],'ninds',[]));
    return
end

hilbert_filtered_lfp = get_hilbert_filtered_call_lfp(session_lfp,filtBank);

for trial_k = 1:n_used_trial
    for t_k = 1:nChunks
        current_lfp = hilbert_filtered_lfp(:,t_idx(t_k,:),used_trial_idx{trial_k},:);
        MIstruct{trial_k,t_k} = calculate_cfc(current_lfp,'selectFilt',selectFilt,...
            'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
        if n_trial_boot > 0
            nChannel = size(hilbert_filtered_lfp,1);
            MIBoot = nan(n_trial_boot,nChannel);
            parfor boot_k = 1:n_trial_boot
                trial_shuffled_lfp = get_trial_shuffled_lfp(current_lfp);
                MI_struct_shuffle = calculate_cfc(trial_shuffled_lfp,'selectFilt',selectFilt,...
                    'n_phase_bins',nbins,'nSurr',nSurr,'randtype',randtype);
                MIBoot(boot_k,:) = [MI_struct_shuffle.MI];
            end
            MIBoot = num2cell(MIBoot,1);
            [MIstruct{trial_k,t_k}.MI_trial_boot] = deal(MIBoot{:});
        end
    end
end


end

function hilbert_filtered_lfp = get_hilbert_filtered_call_lfp(session_lfp,filtBank)

nBat = length(session_lfp);
nFilt = length(filtBank);

nChannel = size(session_lfp{1},1);
nT = size(session_lfp{1},2);
nTrial = size(session_lfp{1},3);

hilbert_filtered_lfp = nan(nChannel,nT,nTrial,nFilt);

for trial_k = 1:nTrial
    for filt_k = 1:nFilt
        if nBat == 1
            current_lfp = squeeze(session_lfp{1}(:,:,trial_k));
        else
            current_lfp = squeeze(session_lfp{filt_k}(:,:,trial_k));
        end
        if ~any(isnan(current_lfp),'all')
            filtered_lfp = filtfilt(filtBank{filt_k},current_lfp');
            hilbert_filtered_lfp(:,:,trial_k,filt_k) = hilbert(filtered_lfp)';
        end
    end
end

end


function trial_shuffled_lfp = get_trial_shuffled_lfp(current_lfp)

trial_shuffled_lfp = current_lfp;
nTrial = size(current_lfp,3);
nFilt = size(current_lfp,4);
for filt_k = nFilt
    permIdx = randperm(nTrial);
    trial_shuffled_lfp(:,:,:,filt_k) = trial_shuffled_lfp(:,:,permIdx,filt_k);
end
end