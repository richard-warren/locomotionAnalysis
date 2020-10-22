function concatTopBotVids(session, varargin)

% old sessions had separate videos for the top and bottom videos of the
% mouse running, but new sessions save a single 'run' video containing both
% views // this repo has been modified to expect the latter // if run, does
% not exist, concatTopBotVids will call ffmpeg to create 'run' video by
% vertically concatenating top and bot vids

% settings
s.bitRate = 10;         % megabits per second
s.moveOriginalTo = '';  % if a directory is provided, a folder named with session name is created there and original runTop and runBot and moved there

% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end
top = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4');
bot = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4');
concat = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4');

fprintf('%s: concatenating runTop.mp4 and runBot.mp4 and saving run.mp4 to disk...\n', session);
command = sprintf('ffmpeg -y -loglevel panic -stats -r 250 -i %s -i %s -filter_complex ''vstack'' -vcodec mpeg4 -vb %iM %s -y', ...
                  top, bot, s.bitRate, concat); %bitRate% ' Evaluation-Tools\Videos\runTopBot.avi';];
[~, ~] = system(command);



% if moveOriginalTo is provided, copy unconcatenated videos to a new directory for backup purposes
if ~isempty(s.moveOriginalTo)
    
    % check that the old and new videos have the same number of frames
    vidTop = VideoReader(top);
    vidCat = VideoReader(concat);
    if vidTop.NumberOfFrames ~= vidCat.NumberOfFrames
        fprintf('FUCK! The concatenated version of %s has a different number of frames!\n', session);
        keyboard
    end
    
    try
        mkdir(fullfile(s.moveOriginalTo, session))
        movefile(top, fullfile(s.moveOriginalTo, session, 'runTop.mp4'))
        movefile(bot, fullfile(s.moveOriginalTo, session, 'runBot.mp4'))
    catch
        fprintf('FUCK! failed to move %s files to new directory!\n', session);
        keyboard
    end
end
disp('all done!')