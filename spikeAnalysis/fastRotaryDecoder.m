file = 'C:\Users\LindseyBuckingham\Google Drive\columbia\analysis\run.mat';
load(file);

%%
aTimes = whEncodA.times';
bTimes = whEncodB.times';
aStates = logical(whEncodA.level)';
bStates = logical(whEncodB.level)';

% sort encoder events by time of occurence
times = [aTimes bTimes];
[times, sortInds] = sort(times);

% create identities vector, where each element represents whether encoder A (0) or B (1) has changed state
% this is the most economical representation of the encoder data, because the states can be inferred at all times if you known the starting state of both channels and the times at which those channels change
identities = [zeros(1, length(aTimes)) ones(1, length(bTimes))];
identities = identities(sortInds);

%%

tic




close all; figure; plot(times, cumsum(deltas)); pimpFig


toc / 60
