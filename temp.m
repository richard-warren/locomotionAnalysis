% PERFORM OPERATIONS ON ALL EPHYS RECORDINGS

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');

for i = 22:height(ephysInfo)
    try
%         plotQualityMetrics(ephysInfo.session{i})
        tic; packContFiles(ephysInfo.session{i}); fprintf('%s: finished in %.1f minutes\n', ephysInfo.session{i}, toc/60)
        showChannelsOverTime(ephysInfo.session{i}, 8)
    catch
        fprintf('problem with session %s\n', ephysInfo.session{i})
    end
end