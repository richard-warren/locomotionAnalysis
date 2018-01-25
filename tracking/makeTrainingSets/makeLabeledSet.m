function makeLabeledSet(className, labeledDataFile, vidFile, subFrameSize, obsPixPositions, posEgs, negEgsPerEg,...
    featureSetting, paws, jitterPixels, jitterNum, maxOverlap, minBrightness)

% !!! need to document, but generally takes hand labeled paw locations and creates features for classifier training
% this is done by taking subframes at paw locations of subFrameSize, and extracting features with getSubFrameFeatures
% it also automatically grabs negative examples in same frames by randomly taking subframes that don't overlap too much with positive examples
% can also create more positive examples by jittering paw locations



% initializations
dataDir = [getenv('OBSDATADIR') 'svm\trainingData\'];
load(labeledDataFile, 'locations', 'locationFrameInds');
nanInds = isnan(locationFrameInds);
locations = locations(:,~nanInds,:);
locationFrameInds = locationFrameInds(~nanInds);
centPad = floor(subFrameSize / 2); % y,x
posEgsCount = 0;
jitterDirections = [1 1; 1 0; 1 -1; 0 1; 0 -1; -1 1; -1 0; -1 -1] * jitterPixels;
numClasses = size(paws,1) + 1; % add one for the not paw class // otherwise, every row in the paws matrix is a class

% sort chronologically (this may make reading video frames faster)
[locationFrameInds, sortInds] = sort(locationFrameInds);
locations = locations(:, sortInds, :);

egsPerFrame = size(locations,3);
imNumberInd = 1;
featureLength = length(getSubFrameFeatures(zeros(subFrameSize), [0 0], featureSetting));

% load video and sample frame
vid = VideoReader(vidFile);
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);


% iterate through frames of all examples (locations)

totalEgs = size(locations,2) * negEgsPerEg * jitterNum;
features = nan(featureLength, totalEgs);
images = nan(prod(subFrameSize), totalEgs); % stores all images in a matrix, so subframes can be viewed prior to feature extraction
labels = nan(1, totalEgs);

for i = randperm(length(locations))
    
    
    % get frame
    frame = rgb2gray(read(vid, locationFrameInds(i)));
    frame = frame - bg;
    
    
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
            
            % get positive examples
            xy = round(locations(1:2, i, j));
            img = getSubFrame(frame, flipud(xy), subFrameSize); % get subframe
            [features(:, imNumberInd), images(:, imNumberInd)] = getSubFrameFeatures(img, xy, featureSetting);

            if ismember(j, paws)
                [row, ~] = ind2sub(size(paws), find(paws==j));
                labels(imNumberInd) = row;
                posEgsCount = posEgsCount+1;
                fprintf('positive eg #%i\n', posEgsCount);
            else
                labels(imNumberInd) = numClasses;
            end

            imNumberInd = imNumberInd+1;



            % get jitered positive examples
            offsetInds = randperm(8);
            offsetInds = offsetInds(1:jitterNum);

            for k = 1:jitterNum

                xyJittered = xy + jitterDirections(offsetInds(k), :)';
                img = getSubFrame(frame, flipud(xyJittered), subFrameSize); % get subframe
                [features(:, imNumberInd), images(:, imNumberInd)] = getSubFrameFeatures(img, xyJittered, featureSetting);

                if ismember(j, paws)
                    [row, ~] = ind2sub(size(paws), find(paws==j));
                    labels(imNumberInd) = row;
                    posEgsCount = posEgsCount+1;
                    fprintf('positive eg #%i\n', posEgsCount);
                else
                    labels(imNumberInd) = numClasses;
                end

                imNumberInd = imNumberInd+1;
            end
                
                
                

            % get negative examples
            for k = 1:negEgsPerEg

                % find a frame that doesn't overlap with positive examples
                acceptableImage = false;

                while ~acceptableImage 

                    pos = [randi([centPad(1)+1 size(frame,1)-centPad(1)-1])...
                           randi([centPad(2)+1 size(frame,2)-centPad(2)-1])]; % y,x
                    temp = egsMask(pos(1)-centPad(1):pos(1)+centPad(1)-1, pos(2)-centPad(2):pos(2)+centPad(2)-1);
                    pixelsOverlap = sum(temp(:));
                    img = getSubFrame(frame, pos, subFrameSize);
                    
                    if (pixelsOverlap/featureLength) < maxOverlap &&...
                       mean(img(:)) > (mean(frame(:))*minBrightness)
                        acceptableImage = true;
                    end
                end

                % store negative example
                [features(:, imNumberInd), images(:, imNumberInd)] = getSubFrameFeatures(img, fliplr(pos), featureSetting);
                labels(imNumberInd) = numClasses;
                imNumberInd = imNumberInd+1;
            end
        end
    end
end

% remove nan values
validInds = ~isnan(labels);
features = features(:,validInds);
images = images(:,validInds);
labels = labels(validInds);


save([dataDir className '\labeledFeatures.mat'], 'features', 'images', 'labels', 'subFrameSize')







