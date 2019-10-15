function spikeRate = getFiringRateDoubleExp(spkTimes, fs, rise, fall, timeLims)
    
    % create kernel
    kernel = arrayfun(@(x) exp(-x/(fall*fs))-exp(-x/(rise*fs)), 0:1*fs);
    kernel = kernel/sum(kernel);
    
    if exist('timeLims', 'var')
        times = timeLims(1) : (1/fs) : timeLims(2);
    else
        times = 0 : (1/fs) : (max(spkTimes)+10*sig);
    end
    
    % convert spike times to binary vector
    spikesVec = hist(spkTimes, times+(1/fs)/2);
    spikesVec = spikesVec.*fs;
    
    % convolve
    spikeRate = conv(spikesVec, kernel, 'same');
end