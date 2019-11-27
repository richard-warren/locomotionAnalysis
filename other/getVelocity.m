function velocity = getVelocity(positions, windowSize, fs)

    % gets velocity at each point in positions by subtracting points to the
    % left and right and then dividing by delta time
    %
    % input         data: (m) position data from which to compute velocity
    %               windowSize: (s) size of window separating samples that are subtracted from one another to compute velocity
    %               fs: (hz) frequency with which data are sampled
    %
    % output        velocity: (m/s) computed velocity
    %
    % note: assumes positions are sampled at even intervals
    

    % compute window size in samples
    windowSmps = round(windowSize*fs);
    if mod(windowSmps,2)==0; windowSmps = windowSmps+1; end % ensure window size is odd

    % compute velocity
    positions = positions(:)';  % force horizontal orientation
    kernel = [1, zeros(1,windowSmps-2), -1];
    velocity = conv(positions, kernel, 'valid') / ((windowSmps-1)/fs);
    pad = (windowSmps-1) / 2;
    velocity = [repelem(velocity(1),pad), velocity, repelem(velocity(end),pad)];
end