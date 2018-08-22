function velocity = getVelocity(data, windowSize, fs)

    % gets velocity by smoothing with first order Savitzky-Golay filter, then differentiating
    %
    % input         data: position data from which to compute velocity
    %               windowSize: size of window separating samples that are subtracted from one another to compute velocity (s)
    %               fs: frequency at which data are sampled
    %
    % output        velocity: computed velocity
    

    % compute window size in samples
    windowSmps = round(windowSize*fs);
    if mod(windowSmps,2)==0; windowSmps = windowSmps+1; end % ensure window size is odd

    % compute velocity
    smoothed = sgolayfilt(data, 1, windowSmps);
    velocity = diff(smoothed) * fs;
    
    % interpolate velocity s.t. values correspond to original times of data
    inds = 1:length(data);
    indsVel = inds(1:end-1) + .5; % values of velocity actually represent velocity in between samples
    velocity = interp1(indsVel, velocity, inds, 'linear', 'extrap');
    

end