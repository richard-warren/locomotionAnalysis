function [positions, times] = rotaryDecoder(aTimes, aStates, bTimes, bStates, encoderSteps, wheelRad, targetFs, session)

% converts quadrature input from rotary encoder into real world positional units (m)
% uses a lookup table method, where every encoder input is represented
% as a transition of a 2 bit state to another 2 bit state (the two bits
% being the states of the two encoder channels). this 4 bit code is
% used to index into a lookup table that tells whether the encoder has
% moved forward, backward, or remained stationary
%
% input        aTimes, bTimes:   times of state shifts in the A and B channels from the rotary encoder
% 	           aStates, bStates: levels of A and B rotary encoder at times of state shifts
%              encoderSteps:     number of steps per rotation in rotary encoder
%              wheelRad:         radius or wheel, or timing pulley, attached to encoder
%              targetFs:         sampling frequency with which data are interpolated (see interpData)
%
% output       positions:        position of wheel in meters
%              times:            times of all position values


% convert inputs to row vectors
aTimes = aTimes(:)';
bTimes = bTimes(:)';
aStates = logical(aStates(:))';
bStates = logical(bStates(:))';

% sort encoder events by time of occurence
times = [aTimes bTimes];
states = [aStates bStates];
[times, sortInds] = sort(times);
states = states(sortInds);

% create 2xn logical matrix, where every column is the state of encoder A at state transitions
allStates = nan(2, length(states));

% this matrix encodes whether a column corresponds to a transition of A (1) or B (2)
identities = [ones(1, length(aTimes)) ones(1, length(bTimes))*2];
identities = identities(sortInds);

% update the changed channel at all state transitions
allStates(1, identities==1) = states(identities==1);
allStates(2, identities==2) = states(identities==2);

% remaining nan values correspond to unchanged channels, so they should
% take the values from the previous state
allStates = fillmissing(allStates, 'previous', 2);

% fill in the final nan in the first column as the opposite of it's subsequent state
nanBins = isnan(allStates(:,1));
allStates(nanBins, 1) = ~allStates(nanBins, 2);

% create transitions matrix, where each column as [prevStateA; prevStateB; currentStateA; currentStateB]
transitions = cat(1, allStates(:,1:end-1), allStates(:,2:end));

% convert columbs to decimal indices and index into lookUp, which converts
% state transitions into deltas
inds = bi2de(transitions')+1;
lookUp = [0,-1,1,0,1,0,0,-1,-1,0,0,1,0,1,-1,0];
deltas = lookUp(inds);
times = times(2:end);

% convert to real-world units (m)
mmPerTic = (2*wheelRad*pi) / encoderSteps;
positions = cumsum(deltas) * (mmPerTic / 1000);


% report whether events are recorded out of order
faultyEventCount = sum(diff(aTimes)<=0) + sum(diff(bTimes)<=0) + length(intersect(aTimes, bTimes));
if faultyEventCount>0
    fprintf('  %s: WARNING: detected %i events simultaneously or out of order!\n', session, faultyEventCount)
    
    % remove duplicate times (necessary for interpolation)
    [times, uniqueInds] = unique(times);
    positions = positions(uniqueInds);
end


% interpolate
[positions, times] = interpData(times, positions, targetFs);








