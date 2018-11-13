

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');

for i = 1:height(ephysInfo)
    disp(ephysInfo.session{i})
    try
        showChannelsOverTime(ephysInfo.session{i}, 8)
    catch
        fprintf('problem with session %s\n', ephysInfo.session{i})
    end
end