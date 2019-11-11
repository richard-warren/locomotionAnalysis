

load(fullfile(getenv('OBSDATADIR'), 'sessions', prev{i}, 'runAnalyzed.mat'), 'obsHeightsVid');


prev = {'190919_000', '190919_001', '190919_002'};
new = {'190921_001', '190921_002', '190921_003'};

prevHgts = cell(1,3);
newHgts = cell(1,3);

for i = 1:3
    load(fullfile(getenv('OBSDATADIR'), 'sessions', prev{i}, 'runAnalyzed.mat'), 'obsHeightsVid');
    prevHgts{i} = obsHeightsVid;
    
    load(fullfile(getenv('OBSDATADIR'), 'sessions', new{i}, 'runAnalyzed.mat'), 'obsHeightsVid');
    newHgts{i} = obsHeightsVid;
end

prevHgts = cat(2, prevHgts{:});
newHgts = cat(2, newHgts{:});

%% check head height across days


% sessions = {'191106_001', '191107_001', '191108_001', '191109_001'};  % sen13
% sessions = {'191106_002', '191107_002', '191108_002', '191109_002'};  % sen14
% sessions = {'191106_003', '191107_003', '191108_003', '191109_000'};  % sen15
sessions = {'191107_004', '191108_000', '191109_003'};  % sen16

for i = 1:length(sessions)
    showWiskHeights(sessions{i});
end







