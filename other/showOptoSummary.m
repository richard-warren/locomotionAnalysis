function showOptoSummary()

% displays a table showing what optogenetics experiments have been done thus far

sessionInfoTemp = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'optoNotes');
sessionInfoTemp = sessionInfoTemp(~cellfun(@isempty, sessionInfoTemp.session) & [sessionInfoTemp.include],:); % remove empty rows, not included sessions
mice = unique(sessionInfoTemp.mouse);
% brainRegions = unique(sessionInfoTemp.brainRegion);
brainRegions = {'alm', 'mtc', 'olf', 'vermis', 'cerInt', 'cerLat'};
summary = cell2table(cell(length(mice), length(brainRegions)+1), 'VariableNames', [{'mouse'} brainRegions]);

for i = 1:length(mice)
    for j = 1:length(brainRegions)
        inds = find(strcmp(sessionInfoTemp.mouse, mice{i}) & strcmp(sessionInfoTemp.brainRegion, brainRegions{j}))';
        entry = '';
        for k = inds; entry = [entry num2str(sessionInfoTemp.mWpeak(k)) upper(sessionInfoTemp.side{k}(1)) ' ']; end
        summary.(brainRegions{j}){i} = entry(1:end-1);
    end
end
summary.mouse = mice;
disp(summary)
