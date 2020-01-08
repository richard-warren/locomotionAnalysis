function pairs = propensityMatching(X, isManip, opts)

% performs propensity score matching // computes propensity score as likelihood that each trial
% came from control or manipulated condition using logistic regression on
% predictors, and finds populations of control and manip trials
% that are matched for propensity score // 'X' is samples x predictors matrix, 
% and isManip is logical vector encoding whether a trial is in
% manipulated condition // pairs is (num matched samples) x (2) matrix, where each row is a matched
% control, manip ind (inds are wrt X and is Manip)

% settings
s.percentileThresh = 10;  % only take trials pairs with propensity score differences in the top percentileThresh percentile
s.predictorNames = {};  % names of columns of X // used to print out how the means of these variables are affected by matching
s.verbose = true;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% fit logit model
glm = fitglm(X, isManip, 'Distribution', 'binomial');
p = predict(glm, X);
ctlInds = find(~isManip);
manipInds = find(isManip);
distances = abs(repmat(p(isManip),1,length(ctlInds)) - p(~isManip)');  % rows are manip trials, cols are control trials, and entries are propensity score differences between each pair of manip/control trials


% greedy search for matched pairs (would be better to replace with hungarian algorithm...)
distancesTemp = distances;
pairNum = floor(s.percentileThresh/100*length(manipInds));  % find the most closely matched pairNum pairs
pairs = nan(pairNum, 2); % numPairs X 2, where columns are control and manip inds

for i = 1:pairNum
    [~, minCol] = min(min(distancesTemp,[],1));
    [~, minRow] = min(distancesTemp(:, minCol));
    pairs(i,:) = [ctlInds(minCol), manipInds(minRow)]; % control ind, manip ind
    distancesTemp(minRow,:) = nan;
    distancesTemp(:,minCol) = nan;
end


% print matched sample predictor means
if s.verbose
    fprintf('\nsamples:  %i/%i', numel(pairs), length(isManip))
    fprintf('\ncontrol:  ')
    fprintf('%.4f  ', nanmean(X(pairs(:,1),:),1));
    fprintf('\nmanip:    ')
    fprintf('%.4f  ', nanmean(X(pairs(:,2),:),1));
    if ~isempty(s.predictorNames); fprintf('\n          '); fprintf('%s ', s.predictorNames{:}); end
    fprintf('\n')
end