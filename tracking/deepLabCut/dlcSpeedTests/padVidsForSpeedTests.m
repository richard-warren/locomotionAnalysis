% warps a video to different dimensions to for mathis DLC benchmarking paper
% goal is to see whether speed varies with squareness of image

% settings
vidName = 'D:\dlcBenchmarking\padtests\original.avi';
heights = 60:4:100;
width = 80;
reps = 10; % reps of each frame (a hack to make videos longer)


vid = VideoReader(vidName);
[folder, name, ext] = fileparts(vidName);

for i = 1:length(heights)
    
    vidWriter = VideoWriter(fullfile(folder, [name sprintf('%ix%ireps%i', heights(i), width, reps) ext]), 'MPEG-4');

    open(vidWriter)
    for j = 1:vid.NumberOfFrames
        
        frame = read(vid,j);
        frame = imresize(frame, [heights(i) width]);
        writeVideo(vidWriter, repmat(frame,1,1,1,reps));
        
    end
    close(vidWriter)
    
end

disp('all done!')

%% scales image to different dimensions

% settings
vidName = 'D:\dlcBenchmarking\sizevsspeed\odor.avi';
dimsNum = 10;
smallestDims = [80 60];
largestDims= [640 480];




% initializations
dims = [round(linspace(smallestDims(1), largestDims(1), dimsNum)); ...
        round(linspace(smallestDims(2), largestDims(2), dimsNum))];


dims = dims + mod(dims,2); % make everything even numbers
[folder, name, ext] = fileparts(vidName);

for i = 1:size(dims,2)
    command = sprintf('ffmpeg -i %s -vb 10M -s %ix%i %s\\%s%04dx%04d.avi', ...
        vidName, dims(1,i), dims(2,i), folder, name, dims(1,i), dims(2,i)); 
    system(command);
end
disp('all done!')











