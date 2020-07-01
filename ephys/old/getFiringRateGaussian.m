function [spikeRate, times] = getFiringRateGaussian(spkTimes, fs, sig, timeLims)

    % create kernel
    kernel = arrayfun(@(x) (1/(sig*sqrt(2*pi))) * exp(-.5*(x/sig)^2), -sig*5:1/fs:sig*5);
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