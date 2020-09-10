% function makeTrialPairVid(sessions, varargin)
% make videos of trials matched for specified vars (all computed via
% getExperimentData) // 'sessions' are the two sessions to match

% temp
sessions = {'180717_005', '180716_003'};
conditions = {'muscimol', 'saline'};

% settings
s.conditions = {};                          % name of conditions for the top and bottom views
s.matchVars = {'trialVel', 'modPawX'};      % vars to match
s.matchConstraints = {'firstModPaw'};       % list of categorical vars that must be equal to eachother in matched trials
s.trials = 10;                              % how many trials to match

% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin

% compute vars for all sessions
data = cell(1,2);
for i = 1:length(sessions)
    sessionData = getExperimentData(sessions{i}, [{'session', 'mouse', 'trial'} s.matchVars s.matchConstraints]);
    data{i} = struct2table(sessionData.data.sessions.trials);
    data{i}.session = repmat(sessions(i), height(data{i}), 1);
    data{i}.condition = ones(height(data{i}),1)*i;
end
data = cat(1, data{:});

% get rid of trials with NaNs
X = table2array(data(:, ismember(data.Properties.VariableNames, [s.matchVars s.matchConstraints])));
data = data(~any(isnan(X),2),:);

%%

X = table2array(data(:, ismember(data.Properties.VariableNames, s.matchVars)));
distances = pdist2(X, X, 'mahalanobis');
% incorporate hard contstraints
for i = 1:length(s.matchConstraints)
    bins = data.(s.matchConstraints{1}) ~= data.(s.matchConstraints{1})';
    distances(bins) = max(distances(:));
end
distances = distances(data.condition==1, data.condition==2);  % maha distance between all condition 1 trials (rows) and all condition 2 trials (cols)




matchedTrials = matchpairs(distances, prctile(distances(:), 90));
inds = sub2ind(size(distances), matchedTrials(:,1), matchedTrials(:,2));  % linear inds for matched trial distances
matchDistances = distances(inds);  % maha distances for each pair of trials
[~, matchInds] = sort(matchDistances);
matchInds = matchInds(1:s.trials);
matchedTrials = matchedTrials(matchInds,:);


% convert back to original inds
c1Inds = find(data.condition==1);
c2Inds = find(data.condition==2);
matchedTrials(:,1) = c1Inds(matchedTrials(:,1));
matchedTrials(:,2) = c2Inds(matchedTrials(:,2));


%% make vids

sessions = cell(2, s.trials);
trials = nan(2, s.trials);
trialText = cell(2, s.trials);
for i = 1:s.trials
    sessions(:,i) = data.session(matchedTrials(i,:));
    trials(:,i) = data.trial(matchedTrials(i,:));
    trialText{1,i} = strjoin(string(table2array(data(matchedTrials(i,1), [s.matchVars s.matchConstraints]))));
    trialText{2,i} = strjoin(string(table2array(data(matchedTrials(i,1), [s.matchVars s.matchConstraints]))));
end

makeTrialPairsVidOld('Z:\loco\obstacleData\papers\hurdles_paper1\movies\matched', sessions, trials, trialText)
















