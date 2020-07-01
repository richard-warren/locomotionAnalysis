function packContFiles(sessions, varargin)

% wrapper function that calls python pack_2 script

% settings
s.pythonPath = 'C:\Users\rick\Anaconda3\envs\deepLabCut\python.exe';
s.highPassFreq = 0; % 0 to skip highpass
s.referencing = 'med';


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if isstr(sessions); sessions = {sessions}; end
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))


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
    ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
    mapFile = ephysInfo.map{strcmp(sessions{i}, ephysInfo.session)};
    load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), 'connected')
%     connected = true(64,1);  % !!! temp (this line sets all channels to be connected)
    connected = [num2str(connected)]'; % string containing binary vector

    % run pack_2
    fprintf('%s: running pack_2...\n', sessions{i})
    fileName = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, ephysFolder);
    commandStr = [s.pythonPath ' ephys\packContFiles.py ' ...,
        fileName ' ' fileNameBase ' ' num2str(fs) ' ' num2str(s.highPassFreq) ' ' s.referencing ' ' connected];
    tic; [~,~] = system(commandStr);
    fprintf('%s: created .dat file in %.1f minutes\n', sessions{1}, toc/60)
end



