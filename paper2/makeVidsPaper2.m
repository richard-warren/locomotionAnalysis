%% make vids sorted by blip amp

load('Z:\loco\obstacleData\data_transfer\to_remote\blip_data.mat', 'blip_data')
[~, sortInds] = sort(blip_data.blip_mean, 'descend');
blip_data = blip_data(sortInds,:);


close all
numSessions = 5;
numTrials = 5;
window = [-.5 1];  % (s) time before and after reward to show

for i = 1:numSessions
    session = blip_data{i, 'session'}{1};
    unit = blip_data{i, 'unit'};
    blipiness = blip_data{i, 'blipiness'}{1};
    [~, sortInds] = sort(blipiness, 'ascend');
    
    inds  = {sortInds(1:numTrials), sortInds(end-numTrials+1:end)};
    names = {'least blippy', 'most blippy'};
    
    for j = 1:2
        windows = blip_data{i, 'rewardTimes'}{1}(sort(inds{j})) + window;
        makeUnitVid(session, unit, ...
            sprintf('Z:\\loco\\obstacleData\\editedVid\\blipSorted\\%s unit %i (%s).avi', session, unit, names{j}), ...
            'vidType', 'showSpecificTimeWindows', 'specificTimeWindows', windows)
    end
    
    
end


%% make hand selected reward vids
session = '201228_000';
unit = 7;

makeUnitVid(session, unit, ...
    sprintf('Z:\\loco\\obstacleData\\editedVid\\blipSorted\\%s unit %i.avi', session, unit), ...
    'vidType', 'showRewardEvents', 'specificRewardTrials', 10:20, 'timeBuffer', [.5 1]);

%% make obstacle vids for most step tuned units

% settings
nunits = 40;
ntrials = 5;
window = [-1.5 1.5];  % time window relative to whisker contact


% load remotely computed step importance
load('Z:\loco\obstacleData\data_transfer\to_remote\dataWithStepImportance.mat', 'data')  % computed on home machine in obstaclestep_sandbox.mat
data = sortrows(data, 'stepImportance', 'descend');


close all
for i = 1:nunits
    % figure out valid times for unit
    neuralData = load(fullfile(getenv('OBSDATADIR'), 'sessions', data.session{i}, 'neuralData.mat'), ...
        'spkRates', 'unit_ids', 'timeStamps');
    unitBin = neuralData.unit_ids==data.unit(i);
    tsub = neuralData.timeStamps(~isnan(neuralData.spkRates(unitBin,:)));
    
    % select trials
    load(fullfile(getenv('OBSDATADIR'), 'sessions', data.session{i}, 'runAnalyzed.mat'), ...
        'wiskContactTimes')
    wiskContactTimes = wiskContactTimes(wiskContactTimes>tsub(1) & wiskContactTimes<tsub(end));
    t = sort(datasample(wiskContactTimes, ntrials, 'replace', false))';
    
    fname = fullfile('Z:\loco\obstacleData\editedVid\stepTunedUnits\', ...
        sprintf('(%03i) %s unit%i.avi', i, data.session{i}, data.unit(i)));
    makeUnitVid(data.session{i}, data.unit(i), fname, ...
        'vidType', 'showSpecificTimeWindows', 'specificTimeWindows', t + window)
end

























