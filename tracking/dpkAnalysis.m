function dpkAnalysis(session, varargin)

% runs deepposekit analysis

s.pythonPath = 'C:\Users\rick\Anaconda3\envs\deepposekit\python.exe';
s.matlabPath = 'C:\Program Files\MATLAB\R2019a\bin\win64';  % this is needed to fix the bug described here: https://www.mathworks.com/matlabcentral/answers/316233-can-t-run-external-program

s.vid = 'run.mp4';
s.model = 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_StackedDenseNet.h5';
s.skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';
s.output = 'trackedFeaturesRaw.csv';

s.verbose = true;



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin

% todo: check for correct video dimensions
fprintf('%s: analyzing %s with deepposekit... ', session, s.vid)
vid = fullfile(getenv('OBSDATADIR'), 'sessions', session, s.vid);
args = [vid ' ' s.model ' ' s.skeleton ' ' s.output];
pathReset = ['set path=%path:' s.matlabPath ';=% &  '];

% run analysis
tic
if s.verbose
    system([pathReset s.pythonPath ' tracking\deepposekit\analyze_video.py ' args]);
else
    [~,~] = system([pathReset s.pythonPath ' tracking\deepposekit\analyze_video.py ' args]);
end
fprintf('deepposekit analysis finished in %.1f minutes\n', toc/60)
