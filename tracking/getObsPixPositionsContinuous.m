function obsPixPositionsContinuous = getObsPixPositionsContinuous(...
    obsPosToWheelPosMappings, wheelTimes, wheelPositions, frameTimeStamps, obsPixPositions, obsPixPositionsUninterped, obsOnTimes, obsOffTimes)


% use wheel position from rotary encoder to infer 'unraveled' obs
% pix positions FOR EACH TRIAL
%
% obsPixPositionsContinuous has one obsPixPositions per trial,
% where row is an 'individual obstacle'. This can be used to
% 'unheadfix' the mouse later, by subtracting these obsPositions
% from the x values of the paw on a trial by trial basis

wheelPositionsInterp = interp1(wheelTimes, wheelPositions, frameTimeStamps); % get position of wheel for all frames
obsPixPositionsContinuous = repmat(obsPixPositions, length(obsOnTimes), 1);

for i = 1:length(obsOnTimes)
    % interp all values expect those where obs is tracked within the trial
    dontInterpBins =  frameTimeStamps>obsOnTimes(i) & ...
                      frameTimeStamps<obsOffTimes(i) & ...
                      ~isnan(obsPixPositionsUninterped); % use obsPixPositions except when it is nan or it it out of frame // otherwise figure it out based on wheel encoder
    if any(~dontInterpBins)
        obsPixPositionsContinuous(i,~dontInterpBins) = wheelPositionsInterp(~dontInterpBins)*obsPosToWheelPosMappings(i,1) + obsPosToWheelPosMappings(i,2);
    end
end