function [positions, times] = motorDecoder(directions, directionTimes, stepTimes, targetFs)

    % converts step and direction signals to stepper motor into real world positional units (m)
    % input     directions:      vector of step directions at times when the direction changes
    %           directionTimes:  times at which directions change 
    %           stepTimes:       times of all steps, regardless of direction
    %           targetFs:        sampling frequency with which data are interpolated (see interpData)
    %
    % output    positions:       position of motor (m) at all stepTimes
    
    % note: the signals necessary to run this function are no longer being recorded as of 190523_000

    
    % stepper motor characteristics
    microStepping = 16;
    motorSteps = 200;
    timingPulleyRad = 15.2789; % mm
    mmPerStep = timingPulleyRad*2*pi / (motorSteps*microStepping);
    
    % convert all vectors to row orientation
    directions = directions(:)';
    directionTimes = directionTimes(:)';
    stepTimes = stepTimes(:)';

    % add initial step direction
    directions = [~directions(1) directions];
    directionTimes = [0 directionTimes];

    % add the final time to directionTimes (make subsequent for loop more elegant)
    directionTimes = [directionTimes max(stepTimes)+1];
    stepDirections = nan(1,length(stepTimes));
    
    
    % iterate through all directional changes to create stepDirections, a vector where all steps are represented with 1 or -1 depending on step direction

    for i = 1:length(directions)

        % determine direction
        if directions(i)
            direction = 1;
        else
            direction = -1;
        end

        % store direction for all steps occuring before the next direction change
        stepDirections(stepTimes>directionTimes(i) & stepTimes<directionTimes(i+1)) = direction;
    end
    
    
    % compute real world positional units (m)
    positions = cumsum(stepDirections) * mmPerStep / 1000;
    
    % interpolate
    [positions, times] = interpData(stepTimes, positions, targetFs);
    
end





