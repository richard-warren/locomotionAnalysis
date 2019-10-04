

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