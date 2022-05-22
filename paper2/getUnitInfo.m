function unitInfo = getUnitInfo(varargin)
% returns table with every recorded unit getting a row // records:
% - mouse
% - session
% - unit
% - target nucleus (nucleus that we were aiming for)
% - nucleus        (nucleus we ended up in, according to reconstruction)
% - ccf location   (3D mm location in allen brain common coordinate framework)


% settings
s.verbose = false;
s.nucleiOnly = true;  % whether to only include units in the cerebellar nuclei
s.frstats = false;    % whether to compute mean and std for each unit

% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% get mice and sessions
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
ephysInfo = ephysInfo(ephysInfo.include==1, :);
mice = unique(ephysInfo.mouse);
unitInfo = cell(1, length(mice));


for i = 1:length(mice)
    sessions = ephysInfo.session(strcmp(ephysInfo.mouse, mice{i}));
    target_nuclei = ephysInfo.target(strcmp(ephysInfo.mouse, mice{i}));
    
    % load registration
    filename = ['E:\lab_files\paper2\histo\registration\' mice{i} '_registration.mat'];
    isRegistered = exist(filename, 'file');
    
    if isRegistered
        load(filename, 'registration');
        tbl = registration(:, {'session', 'unit', 'nucleus', 'ccfMm'});
        tbl = tbl(ismember(tbl.session, sessions), :);  % make sure we are only including sessions that are set to include=1 in ephysInfo
        unitsPerSession = nan(1, length(sessions));
        for j = 1:length(sessions); unitsPerSession(j) = sum(strcmp(registration.session, sessions{j}));end
        target_nuclei = repelem(target_nuclei, unitsPerSession);
%         fprintf('%i %i', length(unique(tbl.session))); disp(length(sessions))
        unitInfo{i} = cat(2, table(repmat(mice(i), height(tbl), 1), target_nuclei(:), ...
            'VariableNames', {'mouse', 'nucleus_target'}), tbl);  % add mouse name column
    else
        % load neural data for each session
        unit_ids = cell(1, length(sessions));
        for j = 1:length(sessions)
            filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{j} '_neuralData.mat']);
            if exist(filename, 'file')
                neuralData = load(filename);
                unit_ids{j} = neuralData.unit_ids;
            else
                fprintf('WARNING! %s does not exist!\n', filename);
            end
        end
        unitsPerSession = cellfun(@length, unit_ids);
        unit_ids = cat(1, unit_ids{:});
        n = length(unit_ids);
        sessions = repelem(sessions, unitsPerSession);
        target_nuclei = repelem(target_nuclei, unitsPerSession);
        unitInfo{i} = table(repmat(mice(i),n,1), target_nuclei(:), sessions(:), unit_ids, repmat({'unregistered'},n,1), nan(n,3), ...
            'VariableNames', {'mouse', 'nucleus_target', 'session', 'unit', 'nucleus', 'ccfMm'});
    end
end

unitInfo = cat(1, unitInfo{:});

% compute mean and std for each unit
if s.frstats
    prevSes = '';
    unitInfo.mean = nan(height(unitInfo), 1);
    unitInfo.std = nan(height(unitInfo), 1);
    
    for i = 1:height(unitInfo)
        
        % load new session data if necessary
        if ~strcmp(prevSes, unitInfo.session{i})
            e = load(['E:\lab_files\paper2\modelling\neuralData\' unitInfo.session{i} '_neuralData.mat']);
            prevSes = unitInfo.session{i};
        end
        unitBin = ismember(e.unit_ids, unitInfo.unit(i));
        unitInfo.mean(i) = nanmean(e.spkRates(unitBin,:));
        unitInfo.std(i) = nanstd(e.spkRates(unitBin,:));
    end
end



if s.verbose
    fprintf('\n\n----------------------------\n')
    fprintf('UNIT COUNTS\n')
    fprintf('----------------------------\n')
    for label = unique(unitInfo.nucleus)'; printCounts(unitInfo, label{1}); end
    fprintf('----------------------------\n\n')
end

if s.nucleiOnly
    unitInfo = unitInfo(ismember(unitInfo.nucleus, {'fastigial', 'interpositus', 'dentate'}),:);
end

end


function printCounts(data, label)
    n = height(data);
    c = sum(strcmp(data.nucleus, label));
    fprintf('%-14s %3i/%i (%2i%%)\n', [label ':'], c, n, round((c/n)*100))
end
