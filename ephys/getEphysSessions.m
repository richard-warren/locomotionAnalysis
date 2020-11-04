function [sessions, neurons] = getEphysSessions()

% returns all sessions in ephysInfo.xlsx for which include==1

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
sessions = ephysInfo.session(ephysInfo.include==1);

if nargout>1
    neurons = cell(1, length(sessions));
    for i = 1:length(sessions)
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{i} '_neuralData.mat']);
        neuralData = load(filename);
        neurons{i} = neuralData.unit_ids;
    end
end
