


files = dir(fullfile(getenv('OBSDATADIR'), 'editedVid', 'reactionTimeVids', '*.csv'));
names = cellfun(@(x) x(24:end-4), {files.name}, 'UniformOutput', false);

spreadSheets = cell(1,length(files));
for i = 1:length(spreadSheets)
    spreadSheets{i} = readtable(fullfile(files(i).folder, files(i).name));
    if i==1; times = nan(height(spreadSheets{1}), length(spreadSheets)); end % initialize matrix on first iteration
    times(:,i) = table2array(spreadSheets{i}(:,4));
end

times(times==999) = nan; % 999 is entered for trials that couldn't be analyzed
times = times*4;


%%

close all;
figure('Color', 'white', 'menubar', 'none');
bins = 2:4:max(times(:))+2;
for i = 1:3
    histogram(times(:,i), bins); hold on
end
legend(names)
set(gca, 'XTick', bins-2, 'box', 'off', 'ycolor', 'white')
xlabel('reaction time (ms)')
