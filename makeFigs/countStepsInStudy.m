% this ridiculous script counts all of the steps in the study

%% count number of trials, from which number of steps can be inferred


experiments = {'baseline', 'mtc_lesion', 'mtc_muscimol', 'senLesion', 'sensoryDependence', 'whiskerTrim'};  % which sheets to analyze

trials = nan(1, length(experiments));
for i = 1:length(experiments)
    fprintf('loading %s data... ', experiments{i})
    load(fullfile(getenv('OBSDATADIR'), 'matlabData', [experiments{i} '_data.mat']), 'data')
    flat = flattenData(data, {'mouse', 'session', 'trial', 'isTrialAnalyzed'});
    trials(i) = length(flat);
    fprintf('%s data loaded...\n', experiments{i})
end
disp('all done!')

%% count ALL steps, not just those occuring in trials

[sessions, experiments] = getAllExperimentSessions();
steps = nan(1,length(sessions));
for i = 1:length(sessions)
    try
        fprintf('\nloading session #%i (%s, %s)... (%.2f)\n', i, sessions{i}, experiments{i}, i/length(sessions))
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'stanceBins')
        steps(i) = max(cumsum(diff(stanceBins(:,1))==1));
    catch
        fprintf('WARNING: %s failed to analyze...\n', sessions{i});
    end
end





