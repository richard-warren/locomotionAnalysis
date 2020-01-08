%% Getting trial numbers for right front paw touch trials

data = getExperimentData('191007_003', {'isBigStep', 'isTrialSuccess', 'isLightOn', 'isVentralContact', 'numTouchFrames'});
Data = data.data.sessions.trials;
flat = flattenData(Data, {'trial','isBigStep', 'isVentralContact', 'isLightOn', 'numTouchFrames', 'paw'});

trials = [flat.trial];
inds = find([flat.paw] == 3 & [flat.isVentralContact] == 1); % 3 is the right front paw
RFpaw_touchTrials = trials(inds);

