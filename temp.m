%%


sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'lesionNotes');
sessionInfo = sessionInfo(strcmp(sessionInfo.mouse, 'sen4'), :);

data = getExperimentData(sessionInfo, 'all');

dv = getDvMatrix(data, 'isTrialSuccess', struct('name', 'sessionNum', 'levels', 1:25), {'mouse'});


%%
close all;
figure('color', 'white', 'MenuBar', 'none', 'Position', [1936 473 1180 321])
sesPlotRick(dv', {'ylabel', 'success rate'});
set(gca, 'YLim', [0 1]);
line([3.5 3.5], [0 1], 'color', 'red')
line([15.5 15.5], [0 1], 'color', 'red')