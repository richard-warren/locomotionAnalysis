function concatTopBotVids(session, varargin)

% old sessions had separate videos for the top and bottom videos of the
% mouse running, but new sessions save a single 'run' video containing both
% views // this repo has been modified to expect the latter // if run, does
% not exist, concatTopBotVids will call ffmpeg to create 'run' video by
% vertically concatenating top and bot vids

% settings
s.bitRate = 10;  % megabits per second

% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end
top = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4');
bot = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4');
concat = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4');

fprintf('%s: concatenating runTop.mp4 and runBot.mp4 and saving run.mp4 to disk...\n', session);
command = sprintf('ffmpeg -y -loglevel panic -stats -r 250 -i %s -i %s -filter_complex ''vstack'' -vcodec mpeg4 -vb %iM %s -y', ...
                  top, bot, s.bitRate, concat); %bitRate% ' Evaluation-Tools\Videos\runTopBot.avi';];
[~, ~] = system(command)