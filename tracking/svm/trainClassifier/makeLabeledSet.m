function makeLabeledSet(className, labeledDataFile, vidFile, subFrameSize,...
    obsPixPositions, posEgs, negEgsPerEg, includeLocations, paws, threshIntensity, jitterPixels, jitterNum)

% !!! need to document


% settings
dataDir = [getenv('OBSDATADIR') 'svm\trainingData\'];
maxOverlap = .25;



% initializations
load(labeledDataFile, 'locations', 'locationFrameInds');
nanInds = isnan(locationFrameInds);
locations = locations(:,~nanInds,:);
locationFrameInds = locationFrameInds(~nanInds);
centPad = floor(subFrameSize / 2); % y,x
posEgsCount = 0;
jitterDirections = [1 1; 1 0; 1 -1; 0 1; 0 -1; -1 1; -1 0; -1 -1] * jitterPixels;

% sort chronologically (this may make reading video frames faster)
[locationFrameInds, sortInds] = sort(locationFrameInds);
locations = locations(:, sortInds, :);

egsPerFrame = size(locations,3);
imNumberInd = 1;
pixPerSub = prod(subFrameSize);

% load video and sample frame
vid = VideoReader(vidFile);
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);


% iterate through frames of all examples (locations)

totalEgs = size(locations,2) * negEgsPerEg * jitterNum;
features = nan(pixPerSub, totalEgs);
labels = nan(1, totalEgs);

for i = randperm(length(locations))
    
    % get frame
    frame = rgb2gray(read(vid, locationFrameInds(i)));
    frame = frame - bg;
    if exist('threshIntensity', 'var')
        frame(frame>threshIntensity) = threshIntensity; % a hack to limit influence of markers shining in bottom view
    end
    
    % mask obstacle
    if ~isnan(obsPixPositions(locationFrameInds(i)))
        frame = maskObs(frame, obsPixPositions(locationFrameInds(i)));
    end
        
    % create mask of locations of positive examples
    egsMask = zeros(size(frame,1), size(frame,2));

    for j = 1:egsPerFrame
        xy = round(locations(1:2, i, j));
        imgInds = {xy(2)-centPad(1):xy(2)+centPad(1), xy(1)-centPad(2):xy(1)+centPad(2)}; % would be smarter to have binary vector keeping track of whether imgInds are valid (if example is too close to edge), so I don't need to compute imgInds multiple times, etc
        imgInds{1}(imgInds{1}<1)=1; imgInds{1}(imgInds{1}>vid.Height)=vid.Height;
        imgInds{2}(imgInds{2}<1)=1; imgInds{2}(imgInds{2}>vid.Width)=vid.Height;
        egsMask(imgInds{1}, imgInds{2}) = 1;
    end


    % save positive and create negative examples
    for j = 1:egsPerFrame
        
        if posEgsCount < posEgs * (1+jitterNum)
            
            xy = round(locations(1:2, i, j));
            [img, isPadded] = getSubFrame(frame, flipud(xy), subFrameSize); % get subframe

            if ~isPadded % if image falls fully within bounds of frame

                if includeLocations
                    img(end, end-1:end) = xy;
                end
                
                features(:, imNumberInd) = img(:);
                
                if ismember(j, paws)
                    labels(imNumberInd) = 1;
                    posEgsCount = posEgsCount+1;
                    fprintf('positive eg #%i\n', posEgsCount);
                else
                    labels(imNumberInd) = 2;
                end
                
                imNumberInd = imNumberInd+1;
                
                
                
                % get jitered positive examples
                offsetInds = randperm(8);
                offsetInds = offsetInds(1:jitterNum);
                
                for k = 1:jitterNum
                    
                    xyJittered = xy + jitterDirections(offsetInds(k), :);
                    [img, isPadded] = getSubFrame(frame, flipud(xyJittered), subFrameSize); % get subframe
                    
                    if ~isPadded
                        
                        if includeLocations
                            img(end, end-1:end) = xyJittered;
                        end

                        features(:, imNumberInd) = img(:);

                        if ismember(j, paws)
                            labels(imNumberInd) = 1;
                            posEgsCount = posEgsCount+1;
                            fprintf('positive eg #%i\n', posEgsCount);
                        else
                            labels(imNumberInd) = 2;
                        end

                        imNumberInd = imNumberInd+1;
                    end
                end
                
                
                

                % create/save negative examples for every positive example
                for k = 1:negEgsPerEg

                    % find a frame that doesn't overlap with positive examples
                    acceptableImage = false;

                    while ~acceptableImage 

                        pos = [randi([centPad(1)+1 size(frame,1)-centPad(1)-1])...
                               randi([centPad(2)+1 size(frame,2)-centPad(2)-1])]; % y,x
                        temp = egsMask(pos(1)-centPad(1):pos(1)+centPad(1)-1, pos(2)-centPad(2):pos(2)+centPad(2)-1);
                        pixelsOverlap = sum(temp(:));
                        [img, isPadded] = getSubFrame(frame, pos, subFrameSize);
                        
                        if (pixelsOverlap/pixPerSub)<maxOverlap && mean(img(:))>mean(frame(:)) && ~isPadded
                            acceptableImage = true;
                        end
                    end

                    % store negative example
                    
                    if includeLocations
                        img(end, end-1:end) = fliplr(pos);
                    end
                    
                    features(:, imNumberInd) = img(:);
                    labels(imNumberInd) = 2;
                    imNumberInd = imNumberInd+1;
                end    
            end
        end
    end
    
    if posEgsCount==posEgs
        break;
    end
end

% remove nan values
validInds = ~isnan(labels);
features = features(:,validInds,:);
labels = labels(validInds);


save([dataDir className '\labeledFeatures.mat'], 'features', 'labels', 'subFrameSize')







