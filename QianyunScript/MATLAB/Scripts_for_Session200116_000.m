%% Getting one big step trials for left front paw
data = getExperimentData('200116_000', {'isBigStep', 'isTrialSuccess', 'isLightOn', 'isLeading', 'firstModPaw'});
Data = data.data.sessions.trials;
flat = flattenData(Data, {'trial','isBigStep', 'isTrialSuccess', 'isLightOn', 'isLeading', 'paw', 'firstModPaw'});

trials = [flat.trial];
inds = find([flat.paw] == 2 & [flat.isLeading] == 1 & [flat.firstModPaw] == 2 & [flat.isBigStep] == 1); % 2 is the left front paw
LFpaw_leadingBigStepTrials = trials(inds);
LFPawLeading_randomTrials = randsample(LFpaw_leadingBigStepTrials, 5)

inds = find([flat.firstModPaw] == 3 & [flat.paw] == 3 & [flat.isBigStep] == 0); % 3 is the right front paw
RFpaw_leadingTwoStepTrials = trials(inds);
randomTrials = randsample(RFpaw_leadingTwoStepTrials, 10)

%% 
data = getExperimentData('200117_000', {'isBigStep', 'isTrialSuccess', 'isLightOn', 'isLeading', 'firstModPaw', 'numTouchFrames', 'isVentralContact'});
Data = data.data.sessions.trials;
flat = flattenData(Data, {'trial','numTouchFrames', 'isVentralContact', 'isTrialSuccess', 'isLightOn', 'isLeading', 'paw', 'firstModPaw'});

trials = [flat.trial];
inds = find([flat.paw] == 1 & [flat.isVentralContact] == 1 & [flat.numTouchFrames] >= 5); % 2 is the left front paw
LHpaw_ventralTouchTrials = [flat(inds).trial];
LHPawVentralTouch_randomTrials = randsample(LHpaw_ventralTouchTrials, 5)



%% Generate videos
session = '200116_000';
unit_id = 52;
fileName = 'Z:\obstacleData\editedVid\vidsWithNeurons\200116_000\unit_52_LFPawLeadingBigStep2.avi'
makeUnitVid(session, unit_id, fileName, {'specificObsTrials', LFPawLeading_randomTrials});

unit_id = 14;
fileName = 'Z:\obstacleData\editedVid\vidsWithNeurons\200116_000\unit_14_randomObsEvent.avi'
makeUnitVid(session, unit_id, fileName);

unit_id = 17;
fileName = 'Z:\obstacleData\editedVid\vidsWithNeurons\200116_000\unit_17_randomObsEvent.avi'
makeUnitVid(session, unit_id, fileName);

unit_id = 50;
fileName = 'Z:\obstacleData\editedVid\vidsWithNeurons\200116_000\unit_50_randomObsEvent.avi'
makeUnitVid(session, unit_id, fileName);

unit_id = 68;
fileName = 'Z:\obstacleData\editedVid\vidsWithNeurons\200116_000\unit_68_randomObsEvent.avi'
makeUnitVid(session, unit_id, fileName);


session = '200117_000';
unit_id = 47;
fileName = 'Z:\obstacleData\editedVid\vidsWithNeurons\200117_000\unit_47_LHPawVentralTouch.avi'
makeUnitVid(session, unit_id, fileName, {'specificObsTrials', LHPawVentralTouch_randomTrials});






