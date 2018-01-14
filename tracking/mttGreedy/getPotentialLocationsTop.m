function potentialLocationsTop = getPotentialLocationsTop(vid, locationsBot, xLinearMapping,...
    model1, model2, subFrameSize1, subFrameSize2, scoreThresh, frameInds, paws, showTracking)

% !!! need to document


% settings
overlapThresh = .6;
yMin = 55; % all pixels below yMin (at the top of the frame) are set to zero in the filtered frame
circRoiPts = [36 172; 224 122; 386 157];
xMaskWidth = 40;


% initializations
xMaskHalfWidth = floor(xMaskWidth/2);
sampleFrame = rgb2gray(read(vid,1));
totalFrames = vid.NumberOfFrames;
kernel = reshape(model1.Beta, subFrameSize1(1), subFrameSize1(2));
wheelMask = double(getWheelMask(circRoiPts, [vid.Height vid.Width]));
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);


% !!! fix x alignment for bottom view
locationsBot = fixTracking(locationsBot);
locationsBot.x = round(locationsBot.x*xLinearMapping(1) + xLinearMapping(2)) + 8; % !!! + 5 is a hack!!!, only temporary



% prepare figure
if showTracking

%     figure; imagesc(-kernel);
    
    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(3,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
    hold on;
    hold on; scatter1 = scatter(rawAxis, 0, 0, 50, [1 1 1], 'filled');
    hold on; scatter2 = scatter(rawAxis, 0, 0, 150, [1 0 0], 'linewidth', 3);
    
    maskAxis = subaxis(3,1,2, 'spacing', 0, 'margin', 0);
    maskIm = image(sampleFrame, 'parent', maskAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
    
    predictAxis = subaxis(3,1,3, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
end


potentialLocationsTop = struct();

for i = frameInds
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = frame - bg;
    frame = getFeatures(frame);

    
   
    
    
    % filter with svm and apply non-maxima suppression
    frameFiltered = -(conv2(double(frame)/model1.KernelParameters.Scale, kernel, 'same') + model1.Bias);
    frameFiltered(frameFiltered < scoreThresh) = 0;
    frameFiltered(1:yMin,:) = 0;
    frameFiltered = frameFiltered .* wheelMask;
    
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

        classes = classify(model2, frameFeatures);
        isPaw = (uint8(classes)==1);

        
        % store data
        potentialLocationsTop(i).x = x(isPaw);
        potentialLocationsTop(i).y = y(isPaw);
        potentialLocationsTop(i).scores = scores(isPaw);
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
        set(scatter1, 'XData', x, 'YData', y);
        set(scatter2, 'XData', x(isPaw), 'YData', y(isPaw));
        
        % pause to reflcet on the little things...
        pause(.2);
%         keyboard
    end
end

close all




