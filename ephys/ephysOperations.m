% PERFORM OPERATIONS ON ALL EPHYS RECORDINGS

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
% ephysInfo = ephysInfo([strcmp(ephysInfo.map, 'D55F')], :); % uncomment to reanalyze data from specific probe
ephysInfo = ephysInfo([strcmp(ephysInfo.spikesSorted, 'yes')], :); % uncomment to reanalyze data that were already spike sorted

for i = 1:height(ephysInfo)
    try
%         formatEphysData(ephysInfo.session{i});
%         plotQualityMetrics(ephysInfo.session{i})
%         tic; packContFiles(ephysInfo.session{i}); fprintf('%s: finished in %.1f minutes\n', ephysInfo.session{i}, toc/60)
%         showChannelsOverTime(ephysInfo.session{i})
        plotManyPSTHs(ephysInfo.session{i}); close all
    catch
        fprintf('problem with session %s\n', ephysInfo.session{i})
    end
end