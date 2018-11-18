% PERFORM OPERATIONS ON ALL EPHYS RECORDINGS

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
% ephysInfo = ephysInfo([strcmp(ephysInfo.map, 'D55F')], :); % uncomment to reanalyze data from specific probe

for i = 3:height(ephysInfo)
    try
%         plotQualityMetrics(ephysInfo.session{i})
%         tic; packContFiles(ephysInfo.session{i}); fprintf('%s: finished in %.1f minutes\n', ephysInfo.session{i}, toc/60)
        showChannelsOverTime(ephysInfo.session{i})
    catch
        fprintf('problem with session %s\n', ephysInfo.session{i})
    end
end