function makeBrainMovie(fig, vidName, varargin)
% make shameless video of brain spinning around in 3D :)

% settings
s.time = 20;            % (s) time for one rotation of the brain
s.fps = 30;            % frames per second
s.startAngle = -45;     % (degrees)
s.height = 35;         % (degrees)


% inits
disp('making brain movie...')
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vidWriter = VideoWriter(vidName);
set(vidWriter, 'FrameRate', s.fps)
set(vidWriter, 'Quality', 50);
open(vidWriter)
numFrames = round(s.time * s.fps);
angles = linspace(s.startAngle, s.startAngle+360, numFrames+1);
angles = angles(1:end-1);
set(gca, 'CameraViewAngleMode', 'Manual');  % prevent axis from rescaling

for i = 1:numFrames
    view([angles(i) s.height]);
    frame = getframe(fig);
    writeVideo(vidWriter, frame.cdata);
end
close(vidWriter)
disp('all done!')