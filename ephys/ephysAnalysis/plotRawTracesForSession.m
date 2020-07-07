function plotRawTracesForSession(session, varargin)


sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), 'unit_ids');

for i = 1:length(unit_ids)
    plotRawTraces(session, unit_ids(i));    
end

end