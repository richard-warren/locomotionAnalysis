function dlcAnalysis(session, varargin)
% runs OLD deeplabcut analysis based on rick's fork of the old deeplabcut repo: https://github.com/rwarren2163/DeepLabCut

% settings
s.verbose = false;
s.dlcPath = 'D:\github\DeepLabCut';
s.output = 'trackedFeatures_run.csv';
s.pythonPath = 'C:\Users\rick\Anaconda3\envs\deepLabCut\python.exe';


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
sesDir = fullfile(getenv('OBSDATADIR'), 'sessions', session);

% concatenate top and bottom views if necessary
if ~exist(fullfile(sesDir, 'run.mp4'), 'file'); concatTopBotVids(newSessions{1}); end  % old sessions were recorded with separate top and bot views, which need to be concatenated

% move video into DLC directory
copyfile(fullfile(sesDir, 'run.mp4'), fullfile(s.dlcPath, 'Evaluation-Tools', 'Videos', 'run.avi'))


% analyze video
fprintf('%s: analyzing run.mp4 with (old) deeplabcut... ', session)
tic
if s.verbose
    system([s.dlcPath(1:2) ' && cd ' fullfile(s.dlcPath, 'Evaluation-Tools') ' && ' s.pythonPath ' AnalyzeVideos.py']); % first move to correct drive, then move to Evaluation-Tools, the run script
else
    [~,~] = system([s.dlcPath(1:2) ' && cd ' fullfile(s.dlcPath, 'Evaluation-Tools') ' && ' s.pythonPath ' AnalyzeVideos.py']);
end
fprintf('%s: (old) deeplabcut analysis finished in %.1f minutes\n', session, toc/60)


% move and rename output file
movefile(fullfile(s.dlcPath, 'Evaluation-Tools', 'Videos', 'trackedFeaturesRaw.csv'), fullfile(sesDir, s.output))

% delete copy of video
delete(fullfile(s.dlcPath, 'Evaluation-Tools', 'Videos', 'run.avi'))

% save metadata
[~, nameNoExtension, ~] = fileparts(s.output);
analysis = 'deeplabcut (old)';
model = 'deeplabcut (old)';
output = s.output;
save(fullfile(sesDir, [nameNoExtension '_metadata.mat']), 'analysis', 'model', 'output')




