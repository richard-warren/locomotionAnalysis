function dpkAnalysis(session, view, varargin)

% runs deepposekit analysis on either run or whisker videos for session //
% 'view' is either 'run' or 'wisk'

s.pythonPath = 'C:\Users\rick\Anaconda3\envs\deepposekit\python.exe';
s.matlabPath = 'C:\Program Files\MATLAB\R2019a\bin\win64';  % these is needed to fix the bug described here: https://www.mathworks.com/matlabcentral/answers/316233-can-t-run-external-program

s.runVid = 'run_originalDimensions.mp4';
s.runModel = 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_StackedDenseNet.h5';
s.runSkeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';
s.runOutput = 'trackedFeatures_run.csv';

s.wiskVid = 'runWisk.mp4';
s.wiskModel = 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_wisk.h5';
s.wiskSkeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';
s.wiskOutput = 'trackedFeatures_wisk.csv';

s.verbose = true;



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin

% todo: check and correct video dimensions


switch view
    case 'run'
        fprintf('%s: running deepposekit on run video... ', session)
        vid = fullfile(getenv('OBSDATADIR'), 'sessions', session, s.runVid);
        args = [vid ' ' s.runModel ' ' s.runSkeleton ' ' s.runOutput];
    case 'wisk'
        fprintf('%s: running deepposekit on whisker video... ', session)
        vid = fullfile(getenv('OBSDATADIR'), 'sessions', session, s.wiskVid);
        args = [vid ' ' s.wiskModel ' ' s.wiskSkeleton ' ' s.wiskOutput];
end
pathReset = ['set path=%path:' s.matlabPath ';=% &  '];

% run analysis
tic
if s.verbose
    system([pathReset s.pythonPath ' tracking\deepposekit\analyze_video.py ' args]);
else
    [~,~] = system([pathReset s.pythonPath ' tracking\deepposekit\analyze_video.py ' args]);
end
fprintf('%s: deepposekit analysis finished in %.1f minutes\n', session, toc/60)
