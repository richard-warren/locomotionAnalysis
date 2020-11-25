function [sessions, neurons, mice] = getEphysSessions()

% returns all sessions in ephysInfo.xlsx for which include==1

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
ephysInfo = ephysInfo(ephysInfo.include==1, :);
sessions = ephysInfo.session;

if nargout>1
    neurons = cell(1, length(sessions));
    for i = 1:length(sessions)
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{i} '_neuralData.mat']);
        neuralData = load(filename);
        neurons{i} = neuralData.unit_ids;
    end
end

if nargout>2
    mice = unique(ephysInfo.mouse);
end