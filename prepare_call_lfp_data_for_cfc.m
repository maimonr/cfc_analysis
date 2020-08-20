function [session_lfp_ds, call_bat_nums, used_call_IDs, lfp_call_offset] = prepare_call_lfp_data_for_cfc(current_lfp_fname,activeChannels,varargin)

pnames = {'avgTetrode','ds_factor','tRange','artifact_nStd_factor','max_artifact_frac','baselineRange'};
dflts  = {true,2,[],5,0.05,[-Inf -2]};
[avgTetrode,ds_factor,tRange,artifact_nStd_factor,max_artifact_frac,baselineRange] = internal.stats.parseArgs(pnames,dflts,varargin{:});

call_trig_csc = load(current_lfp_fname,'call_trig_csc','batNums','lfp_call_offset','used_call_IDs','active_channels');
call_bat_nums = call_trig_csc.batNums;
used_call_IDs = call_trig_csc.used_call_IDs;
lfp_call_offset = call_trig_csc.lfp_call_offset;
[activeChannels,channelIdx] = intersect(call_trig_csc.active_channels,activeChannels);

call_trig_csc = call_trig_csc.call_trig_csc(:,:,channelIdx);

if avgTetrode
    nChannel = 4;
    nChannel_per_tt = 4;
else
    nChannel = size(call_trig_csc,3);
end

nT = size(call_trig_csc,1);
t = linspace(-lfp_call_offset,lfp_call_offset,nT);
if ~isempty(tRange)
    [~,t_idx] = inRange(t,tRange);
    call_trig_csc = call_trig_csc(t_idx,:,:,:);
    nT = sum(t_idx);
end

nTrial = size(call_trig_csc,2);

if nTrial < 1
    [session_lfp_ds, call_bat_nums] = deal(NaN);
    return
end

if max_artifact_frac > 0
    call_trig_csc = remove_lfp_artifacts(call_trig_csc,t,artifact_nStd_factor,max_artifact_frac,baselineRange);
end

session_lfp = zeros(nChannel,nT,nTrial);

if avgTetrode
    tetrodeChannels = reshape(0:(nChannel*nChannel_per_tt)-1,nChannel,[]);
    for tt_k = 1:nChannel
        used_channel_idx = ismember(activeChannels,tetrodeChannels(:,tt_k));
        current_lfp = squeeze(nanmean(call_trig_csc(:,:,used_channel_idx),3));
        session_lfp(tt_k,:,:) = current_lfp;
    end
else
    session_lfp = permute(call_trig_csc,[3 1 2]);
    if nChannel ~= length(activeChannels) && nChannel == 16
        session_lfp = session_lfp(activeChannels,:,:);
    end
end

session_lfp_ds = nan(nChannel,round(nT/ds_factor),nTrial);
for ch_k = 1:nChannel
    lfp_tmp = squeeze(session_lfp(ch_k,:,:));
    session_lfp_ds(ch_k,:,:) = resample(lfp_tmp,1,ds_factor);
end

end

function call_trig_csc = remove_lfp_artifacts(call_trig_csc,t,artifact_nStd_factor,max_artifact_frac,baselineRange)

nT = size(call_trig_csc,1);
nTrial = size(call_trig_csc,2);
nChannel = size(call_trig_csc,3);

[~,baselineIdx] = inRange(t,baselineRange);
mu = mean(call_trig_csc(baselineIdx,:,:),[1 2]);
sigma = std(call_trig_csc(baselineIdx,:,:),[],[1 2]);

for trial_k = 1:nTrial
    for ch_k = 1:nChannel
        artifactFrac = sum(abs(call_trig_csc(:,trial_k,ch_k) - mu) > artifact_nStd_factor*sigma)/nT;
        if artifactFrac > max_artifact_frac
            call_trig_csc(:,trial_k,ch_k) = NaN;
        end
    end
end

end