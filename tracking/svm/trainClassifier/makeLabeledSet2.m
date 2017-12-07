function makeLabeledSet2(className, labeledDataFile, vidFile)

% !!! need to document


% settings
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trainingImages\';
subFrameSize = [40 40]; % y,x
negEgsPerEg = 5;
maxOverlap = .01;


% initializations
load(labeledDataFile, 'locations', 'locationFrameInds');
nanInds = isnan(locationFrameInds);
locations = locations(:,~nanInds,:);
locationFrameInds = locationFrameInds(~nanInds); % !!! sorting these chhronologically will likely make getting vid frames faster
centPad = floor(subFrameSize / 2); % y,x

egsPerFrame = size(locations,3);
imNumberInd = 1;
negImNumberInd = 1;
pixPerSub = prod(subFrameSize);

% load video and sample frame
vid = VideoReader(vidFile);
bg = getBgImage(vid, 1000, false);


% iterate through all example frame
for i = 1:length(locations)
    
    % get frame
    frame = rgb2gray(read(vid, locationFrameInds(i)));
    frame = frame - bg;
        
    % create mask of locations of positive examples
    egsMask = zeros(size(frame,1), size(frame,2));

    for j=1:egsPerFrame
        xy = round(locations(1:2, i, j));
        imgInds = {xy(2)-centPad(1):xy(2)+centPad(1), xy(1)-centPad(2):xy(1)+centPad(2)}; % would be smarter to have binary vector keeping track of whether imgInds are valid (if example is too close to edge), so I don't need to compute imgInds multiple times, etc
        imgInds{1}(imgInds{1}<1)=1; imgInds{1}(imgInds{1}>vid.Height)=vid.Height;
        imgInds{2}(imgInds{2}<1)=1; imgInds{2}(imgInds{2}>vid.Width)=vid.Height;
        egsMask(imgInds{1}, imgInds{2}) = 1;
    end


    % save positive and create negative examples
    for j=1:egsPerFrame
        
        xy = round(locations(1:2, i, j));
        imgInds = {xy(2)-centPad(1):xy(2)+centPad(1), xy(1)-centPad(2):xy(1)+centPad(2)};
        
        if ~any(imgInds{1}<1 | imgInds{1}>vid.Height) && ~any(imgInds{2}<1 | imgInds{2}>vid.Width)
            img = frame(imgInds{1}, imgInds{2});
            save([dataDir className '\positive\img' num2str(imNumberInd) '.mat'], 'img');
            imNumberInd = imNumberInd+1;
            fprintf('positive eg #%i\n', imNumberInd);

            % create/save negative examples for every positive example
            for k=1:negEgsPerEg

                % find a frame that doesn't overlap with positive examples
                acceptableImage = false;

                while ~acceptableImage 

                    pos = [randi([centPad(1)+1 size(frame,1)-centPad(1)-1])...
                           randi([centPad(2)+1 size(frame,2)-centPad(2)-1])]; % y,x
                    temp = egsMask(pos(1)-centPad(1):pos(1)+centPad(1), pos(2)-centPad(2):pos(2)+centPad(2));
                    pixelsOverlap = sum(temp(:));
                    img = frame(pos(1)-centPad(1):pos(1)+centPad(1), pos(2)-centPad(2):pos(2)+centPad(2));
                    
                    if (pixelsOverlap/pixPerSub)<maxOverlap && mean(img(:))>mean(frame(:))
                        acceptableImage = true;
                    end
                end

                % save negative example
                save([dataDir className '\negative\img' num2str(negImNumberInd) '.mat'], 'img');
                negImNumberInd = negImNumberInd+1;
            end    
        end
    end
end










