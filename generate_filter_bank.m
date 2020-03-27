function filtBank = generate_filter_bank(filtType,band_pass_freqs,stop_band_freqs,attenuationDb,pass_ripple,fs)

nFilts = size(band_pass_freqs,1);
filtBank = cell(1,nFilts);

for filt_k = 1:nFilts
    
    filtBank{filt_k} = designfilt(filtType, 'StopbandFrequency1', stop_band_freqs(filt_k,1), 'PassbandFrequency1', band_pass_freqs(filt_k,1),...
        'PassbandFrequency2', band_pass_freqs(filt_k,2), 'StopbandFrequency2', stop_band_freqs(filt_k,2), 'StopbandAttenuation1', attenuationDb(filt_k,1),...
        'PassbandRipple', pass_ripple(filt_k), 'StopbandAttenuation2', attenuationDb(filt_k,2), 'SampleRate', fs);
end
end