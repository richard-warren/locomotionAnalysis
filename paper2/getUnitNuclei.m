function data = getUnitNuclei()
% returns table with nuclei for each recorded cell in each session // also
% reports whether histo has yet to be analyzed // does this by loading
% mouse_registration.mat files


ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
ephysInfo = ephysInfo(ephysInfo.include==1, :);
mice = unique(ephysInfo.mouse);
data = cell(1,length(mice));

for i = 1:length(mice)
    sessions = ephysInfo.session(strcmp(ephysInfo.mouse, mice{i}));
    
    % load registration
    filename = ['E:\lab_files\paper2\histo\registration\' mice{i} '_registration.mat'];
    isRegistered = exist(filename, 'file');
    
    if isRegistered
        load(filename, 'registration');
        tbl = registration(:, {'session', 'unit', 'nucleus'});
        data{i} = cat(2, table(repmat(mice(i), height(tbl), 1), ...
            'VariableNames', {'mouse'}), tbl);  % add mouse name column
    else
        % load neural data for each session
        unit_ids = cell(1, length(sessions));
        for j = 1:length(sessions)
            neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{j} '_neuralData.mat']));
            unit_ids{j} = neuralData.unit_ids;
        end
        unitsPerSession = cellfun(@length, unit_ids);
        unit_ids = cat(1, unit_ids{:});
        n = length(unit_ids);
        sessions = repelem(sessions, unitsPerSession);
        data{i} = table(repmat(mice(i),n,1), sessions(:), unit_ids, repmat({'unregistered'},n,1), ...
            'VariableNames', {'mouse', 'session', 'unit', 'nucleus'});
    end
end

data = cat(1, data{:});

fprintf('\n\nUNIT COUNTS\n')
fprintf('-----------\n')
for label = unique(data.nucleus)'; printCounts(data, label{1}); end
% fprintf('fastigial:    %3i/%i\n', sum(strcmp(data.nucleus, 'fastigial')), n)
% fprintf('interpositus: %3i/%i\n', sum(strcmp(data.nucleus, 'interpositus')), n)
% fprintf('dentate:      %3i/%i\n', sum(strcmp(data.nucleus, 'dentate')), n)
% fprintf('other:        %3i/%i\n', sum(strcmp(data.nucleus, 'none')), n)
% fprintf('unregistered: %3i/%i\n', sum(strcmp(data.nucleus, 'unregistered')), n)

end


function printCounts(data, label)
    n = height(data);
    c = sum(strcmp(data.nucleus, label));
    fprintf('%-14s %3i/%i (%2i%%)\n', [label ':'], c, n, round((c/n)*100))
end
