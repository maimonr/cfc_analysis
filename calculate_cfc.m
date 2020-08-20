function  MIstruct = calculate_cfc(filteredHilb,varargin)

pnames = {'selectFilt','n_phase_bins','nSurr','randtype'};
dflts  = {[],18,0,2};
[selectFilt,nbins,nSurr,randtype] = internal.stats.parseArgs(pnames,dflts,varargin{:});

miParams = struct('nbins',nbins,'nSurr',nSurr,'randtype',randtype);

hilbAmp = abs(filteredHilb(:,:,:,selectFilt(:,2)));
hilbPhase = angle(filteredHilb(:,:,:,selectFilt(:,1)));
MIstruct_tmp = wrap_get_mi(hilbAmp,hilbPhase,miParams);

MIstruct_tmp = reshape(vertcat(MIstruct_tmp{:}),size(MIstruct_tmp));
MI_struct_fieldnames = fieldnames(MIstruct_tmp);

MIstruct = struct('MI',[],'MIp',[],'MIsurrmi',[],'MIsurr',[],'mean_amps',[],'ninds',[]);

for f_k = 1:length(MI_struct_fieldnames)
    fName = MI_struct_fieldnames{f_k};
    current_data = arrayfun(@(x) squeeze(x.(fName)),MIstruct_tmp,'un',0);
    current_ndims = ndims(current_data{1});
    if sum(size(current_data) ~=1) == 1
        current_data = cellfun(@(x) reshape(x,length(x),1),current_data,'un',0);
        current_ndims = 1;
    end
    MIstruct.(fName) = cat(current_ndims+1,current_data{:});
end

end

function MIstruct = wrap_get_mi(hilbAmp,hilbPhase,miParams)

nChannel = size(hilbAmp,1);
nFilt = size(hilbAmp,4);
MIstruct = cell(nChannel,nFilt);

for filt_k = 1:nFilt
    current_hilb_amp = hilbAmp(:,:,:,filt_k);
    current_hilb_phase = hilbPhase(:,:,:,filt_k);
    for ch_k = 1:nChannel
        miInput_amp = reshape(squeeze(current_hilb_amp(ch_k,:,:)),1,[])';
        miInput_phase = reshape(squeeze(current_hilb_phase(ch_k,:,:)),1,[])';
        MIstruct{ch_k,filt_k} = get_mi(miInput_phase,miInput_amp,miParams.nbins,miParams.nSurr,miParams.randtype);
    end
end


end