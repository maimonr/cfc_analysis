function  MIstruct = calculate_cfc(filteredHilb,varargin)

pnames = {'selectFilt','n_phase_bins','nSurr','randtype'};
dflts  = {[],18,0,2};
[selectFilt,nbins,nSurr,randtype] = internal.stats.parseArgs(pnames,dflts,varargin{:});

miParams = struct('nbins',nbins,'nSurr',nSurr,'randtype',randtype);

hilbAmp = abs(filteredHilb(:,:,:,selectFilt(:,2)));
hilbPhase = angle(filteredHilb(:,:,:,selectFilt(:,1)));
MIstruct = wrap_get_mi(hilbAmp,hilbPhase,miParams);

MIstruct = reshape(vertcat(MIstruct{:}),size(MIstruct));

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