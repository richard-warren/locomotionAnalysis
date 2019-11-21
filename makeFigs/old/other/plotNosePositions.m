sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'sessions');
sessionInfo = sessionInfo(348:end,:);
%%
mice = unique(sessionInfo.mouse);
[nosePosits, wheelCenters] = deal(cell(1,length(mice)));

for i = 1:length(mice)
    
    mouseInds = find(strcmp(sessionInfo.mouse, mice{i}));
    mouseNosePosits = nan(size(mouseInds));
    mouseWheelCenters = nan(2, length(mouseInds));
    fprintf('%s: getting data for session ', mice{i})
    
    for j = 1:length(mouseInds)
        try
            session = sessionInfo.session{mouseInds(j)};
            data = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
                'nosePos', 'wheelCenter');
            mouseNosePosits(j) = data.nosePos(1);
            mouseWheelCenters(:,j) = data.wheelCenter;
            fprintf('%i ', j)
        end
    end
    nosePosits{i} = mouseNosePosits;
    wheelCenters{i} = mouseWheelCenters;
    fprintf('\n')
end



%%
close all; figure;

for i = 1:length(mice)
    histogram(nosePosits{i}, 375:385); hold on
end