

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');

for i = 17:height(ephysInfo)
    try
%         showChannelsOverTime(ephysInfo.session{i}, 8)
        plotQualityMetrics(ephysInfo.session{i})
    catch
        fprintf('problem with session %s\n', ephysInfo.session{i})
    end
end