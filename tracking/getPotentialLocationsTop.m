function potentialLocationsTop = getPotentialLocationsTop(vid, locationsBot,...
    model1, model2, classNum, subFrameSize1, subFrameSize2, scoreThresh, frameInds, paws, showTracking)

% !!! need to document


% settings
overlapThresh = .6;
yMin = 55; % all pixels below yMin (at the top of the frame) are set to zero in the filtered frame
circRoiPts = [36 155; 212 103; 374 132];
xMaskWidth = 40;


% initializations
xMaskHalfWidth = floor(xMaskWidth/2);
sampleFrame = rgb2gray(read(vid,1));
totalFrames = vid.NumberOfFrames;
kernel = reshape(model1.Beta, subFrameSize1(1), subFrameSize1(2));
% wheelMask = double(getWheelMask(circRoiPts, [vid.Height vid.Width]));
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);
cmap = hsv(classNum);




% prepare figure
if showTracking

%     figure; imagesc(-kernel);
    
    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(3,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
    hold on;
    hold on; scatterAll = scatter(rawAxis, 0, 0, 50, [1 1 1], 'filled');
    
    scatterPaws = cell(1, classNum); % shows results of second round of classification
    for i = 1:classNum
        hold on; scatterPaws{i} = scatter(rawAxis, 0, 0, 150, cmap(i,:), 'linewidth', 3);
    end
    
    maskAxis = subaxis(3,1,2, 'spacing', 0, 'margin', 0);
    maskIm = image(sampleFrame, 'parent', maskAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
    
    predictAxis = subaxis(3,1,3, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
end


potentialLocationsTop(max(frameInds)) = struct();

for i = frameInds
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = frame - bg;
    
    
    % filter with svm and apply non-maxima suppression
    frameFiltered = -(conv2(double(frame)/model1.KernelParameters.Scale, kernel, 'same') + model1.Bias);
    
    frameFiltered(frameFiltered < scoreThresh) = scoreThresh;
    frameFiltered = frameFiltered - scoreThresh;
%     frameFiltered(1:yMin,:) = 0;
%     frameFiltered = frameFiltered .* wheelMask;
    
     % mask x positions out of range
%     xMask = double(zeros(size(frame)));
%     
%     for j = paws%1:4
%         if ~isnan(locationsBot.x(i,j))
%             
%             % get mask indices for single paw
%             inds = locationsBot.x(i,j)-xMaskHalfWidth : locationsBot.x(i,j)+xMaskHalfWidth;
%             inds(inds<1) = 1;
%             inds(inds>vid.Width) = vid.Width;
%             
%             % incorporate paw mask into mask
%             xMask(:,inds) = 1;
%         end
%     end
%     frameFiltered = frameFiltered .* xMask;
%     
    [x, y, scores] = nonMaximumSupress(frameFiltered, subFrameSize1, overlapThresh);
    
    if ~isempty(x)
        
        % perform second round of classification (cnn)
        dims = model2.Layers(1).InputSize;
        frameFeatures = nan(dims(1), dims(2), 3, length(x));
        for j = 1:length(x)
            img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
            img = uint8(imresize(img, 'outputsize', model2.Layers(1).InputSize(1:2)));
            img = repmat(img, 1, 1, 3);
            frameFeatures(:,:,:,j) = img;
        end

        classes = uint8(classify(model2, frameFeatures));

        
        % store data
        pawBins = classes<=classNum;
        potentialLocationsTop(i).x = x(pawBins);
        potentialLocationsTop(i).z = y(pawBins);
        potentialLocationsTop(i).scores = scores(pawBins);
        potentialLocationsTop(i).class = classes(pawBins);
    end
    
    
    if showTracking
        
        % put lines in top frame
        for j = paws%1:4
            if locationsBot.x(i,j)>0 && locationsBot.x(i,j)<vid.Width
                frame(:,locationsBot.x(i,j)) = 255;
            end
        end
        
        % update figure
        set(rawIm, 'CData', frame);
%         set(maskIm, 'CData', frameMasked);
        set(predictIm, 'CData', frameFiltered)
        set(scatterAll, 'XData', x, 'YData', y);
        
        for j=1:classNum
            set(scatterPaws{j}, 'XData', x(classes==j), 'YData', y(classes==j));
        end
        
        % pause to reflcet on the little things...
        pause(.2);
    end
end

close all




