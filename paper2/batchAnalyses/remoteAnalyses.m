%% perform analyses that require access to engram, which is slow from home ethernet

sessions = getEphysSessions();
overwrite = false;
dt = .01;

tic
for i = 1:length(sessions)
%     try
        % format ephys data
        filename = fullfile(getenv('OBSDATADIR'), 'data_transfer', 'neuralData', [sessions{i} '_neuralData.mat']);
        if overwrite || ~exist(filename, 'file')
            formatEphysData(sessions{i}, 'outputFileName', filename, 'kernel', 'gauss', 'kernelSig', .02, 'plot', false)
        end
        
        % predictors
        filename = fullfile(getenv('OBSDATADIR'), 'data_transfer', 'predictors', [sessions{i} '_predictors.mat']);
        if overwrite || ~exist(filename, 'file')
            getPredictors(sessions{i}, 'outputFileName', filename, 'plot', true, 'visible', 'off', 'dt', dt)
        end
        
        % experimentData (for single session)
        filename = fullfile(getenv('OBSDATADIR'), 'data_transfer', 'expData', [sessions{i} '_expData.mat']);
        if overwrite || ~exist(filename, 'file')
            getExperimentData(sessions{i}, 'all', 'outputFileName', filename);
        end
        
        % stepData (for single session)
        filename = fullfile(getenv('OBSDATADIR'), 'data_transfer', 'stepData', [sessions{i} '_stepData.mat']);
        if overwrite || ~exist(filename, 'file')
            getStepData(sessions{i}, 'outputFileName', filename);
        end
            
%     catch exception
%         fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
%     end
end 
fprintf('finished in %.1f minutes\n', toc/60)

%% histo

% find out which brains are ready
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'), 'ephysHistoTable')
mice = unique(ephysHistoTable.mouseID);

overwrite = false;

for i = 1:length(mice)
    
    % prepare histo labels
    filename = fullfile(getenv('OBSDATADIR'), 'data_transfer', 'histoLabels', [mice{i} '_histoLabels.mat']);
    if ~exist(filename, 'file') || overwrite
        prepareHistoLabels(mice{i}, filename)
    end
    pause(.1)  % allows figure to render properly
end

%% make paw contact vids

% load (remotely computed) paw contact data
load('Z:\loco\obstacleData\data_transfer\to_remote\pawcontact.mat', 'data')


% settings
% session = '181004_003'; unit = 74;
session = '191007_003'; unit = 20;
% session = '200819_000'; unit = 192;
% session = '200130_000'; unit = 14;

paw = 3;
touch = 'dorsal';  % 'dorsal' or 'ventral'
tlims = [-.4 .4];  % (s) time pre and post touch to show


% get touch times
touches = data{session, [touch '_touches']}{paw};
tlims = tlims + touches;

fprintf('total touches: %i\n', length(touches))
filename = fullfile('C:\Users\rick\Desktop\temp_vids\', sprintf('%sunit%ipaw%i%s.avi', session, unit, paw, touch));
makeUnitVid(session, unit, filename, 'specificTimeWindows', tlims, 'vidType', 'showSpecificTimeWindows')  % make vid with unit
% makeVid(filename, session, 'tlims', tlims, 'includeWiskCam', false)  % make vid with behavior only

































