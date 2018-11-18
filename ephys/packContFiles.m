function packContFiles(sessions)

% wrapper function that calls python pack_2 script

% settings
pythonPath = 'C:\Users\rick\Anaconda3\python.exe';
highPassFreq = 300;
referencing = 'med';

% initializations
if isstr(sessions); sessions = {sessions}; end

for i = 1:length(sessions)
    
    % get session data
    files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}));
    ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;
    contFiles = dir(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, ephysFolder, '*.continuous'));
    contFiles = contFiles(~contains({contFiles.name}, 'AUX')); % remove AUX channels
    fileNameBase = contFiles(1).name(1:3);
    [~, ~, info] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, ephysFolder, [fileNameBase '_CH1.continuous']));
    fs = info.header.sampleRate;
    
    % get channel mapping
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
    ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
    mapFile = ephysInfo.map{strcmp(sessions{i}, ephysInfo.session)};
    load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), 'connected')
    connected = [num2str(connected)]'; % string containing binary vector

    % run pack_2
    fprintf('%s: running pack_2...\n', sessions{i})
    fileName = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, ephysFolder);
    commandStr = [pythonPath, ' ephys\packContFiles.py ' ...,
        fileName ' ' fileNameBase ' ' num2str(fs) ' ' num2str(highPassFreq) ' ' referencing ' ' connected];
    [~,~] = system(commandStr);
end