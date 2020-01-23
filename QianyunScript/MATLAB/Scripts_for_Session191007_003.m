%% Getting trial numbers for right front paw touch trials

data = getExperimentData('191007_003', {'isBigStep', 'isTrialSuccess', 'isLightOn', 'isVentralContact', 'numTouchFrames'});
Data = data.data.sessions.trials;
flat = flattenData(Data, {'trial','isBigStep', 'isVentralContact', 'isLightOn', 'numTouchFrames', 'paw'});

trials = [flat.trial];
inds = find([flat.paw] == 3 & [flat.isVentralContact] == 1); % 3 is the right front paw
RFpaw_touchTrials = trials(inds);

%% Getting one big step trials for left front paw
data = getExperimentData('191007_003', {'isBigStep', 'isTrialSuccess', 'isLightOn', 'isLeading', 'firstModPaw'});
Data = data.data.sessions.trials;
flat = flattenData(Data, {'trial','isBigStep', 'isTrialSuccess', 'isLightOn', 'isLeading', 'paw', 'firstModPaw'});

trials = [flat.trial];
inds = find([flat.paw] == 2 & [flat.isLeading] == 1 & [flat.firstModPaw] == 2 & [flat.isBigStep] == 1); % 2 is the left front paw
LFpaw_leadingBigStepTrials = trials(inds);
randomTrials = randsample(LFpaw_leadingBigStepTrials, 10)

inds = find([flat.firstModPaw] == 3 & [flat.paw] == 3 & [flat.isBigStep] == 0); % 3 is the right front paw
RFpaw_leadingTwoStepTrials = trials(inds);
randomTrials = randsample(RFpaw_leadingTwoStepTrials, 10)
