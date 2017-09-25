function [positions, times] = rotaryDecoder(aTimes, aStates, bTimes, bStates)

    % converts quadrature input from rotary encoder into real world positional units (m)
    %
    % input        aTimes, bTimes:   times of state shifts in the A and B channels from the rotary encoder
    % 	           aStates, bStates: levels of A and B rotary encoder at times of state shifts
    % output       positions:        position of wheel in meters
    %              times:            times of all position values
    

    tic
    % wheel, encoder settings
    encoderSteps = 2880; % 720cpr * 4
    wheelRad = 95.25; % mm
    mmPerTic = (2*wheelRad*pi) / encoderSteps;
    
    % convert all inputs to row vectors
    aTimes = aTimes(:)';
    bTimes = bTimes(:)';
    aStates = logical(aStates(:))';
    bStates = logical(bStates(:))';

    % sort all time indices
    times = [aTimes bTimes];
    [times, timesSortedInds] = sort(times);

    % create 2xN matrix of encoder states, sorted with time indices determind above
    % each row is an encoder (A or B), and NaN values temporarily stand in for times when the other encoder changes state
    allStates = nan(2, length(aTimes)+length(bTimes));

    temp = [aStates, nan(1, length(bStates))];
    allStates(1,:) = temp(timesSortedInds);

    temp = [nan(1, length(aStates)), bStates];
    allStates(2,:) = temp(timesSortedInds);

    
    
    % replace all NaN values with opposite of subsequent real number
    for i = 1:2
        for j = find(isnan(allStates(i,:)))

            % find next non-NaN value and store the opposite of that value in current index
            nextValInd = find(~isnan(allStates(i,j:end)), 1, 'first');
            nextVal = allStates(i, j+nextValInd-1);

            if ~isempty(nextVal)
                allStates(i,j) = ~nextVal;
            else
                % fill in missing NaN values at the very end
                allStates(i,j:end) = allStates(i,j-1);
            end
        end    
    end
    
    

    % convert state transitions into code for forward (1) or backward (-1) tics using lookup table method
    % lookup table method described here: http://makeatronics.blogspot.com/2013/02/efficiently-reading-quadrature-with.html
    
    positions = nan(1, length(times));
    lookUp = [0,-1,1,0,1,0,0,-1,-1,0,0,1,0,1,-1,0]; % multiply by -1 to make positions positive rather than negative

    % compute first position
    % (the computation for the first position is a little different because baseline state of encoder pins must be inferred from level at first state shift)
    transitionCode = [num2str(~aStates(1))...
                      num2str(~bStates(1))... % the initial state is always the opposite of the state detected at the first event
                      num2str(allStates(1,1))...
                      num2str(allStates(2,1))];
    positions(1) = lookUp(bin2dec(transitionCode) + 1); % add 1 to compensate because matlab uses 1 based indexing

    
    % compute remaining positions
    for i = 2:length(times)

        transitionCode = [num2str(allStates(1,i-1))...
                          num2str(allStates(2,i-1))...
                          num2str(allStates(1,i))...
                          num2str(allStates(2,i))];

        positions(i) = positions(i-1) + lookUp(bin2dec(transitionCode) + 1); % add 1 to compensate because matlab uses 1 based indexing
    end
    
%     keyboard

    % convert to real-world units (m)
    positions = positions * (mmPerTic / 1000);
    
    decodingTime = toc/60; % minutes
    fprintf('rotary decoding time: %f minutes\n', decodingTime)
    
end











