function makeLabeledSet(className, imNumbers, egsPerFrame, file)

    % TO DO:
    % - document
    % - save position of egs, make that a component of th features?
    % - make imNumbers not necessarily divisible by egsPerFrame
    % - add box position constraint function
    
    % user settings
    dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trainingImages\';
%     subHgtWid = [40 40]; % use this for pawBot
    subHgtWid = [60 40]; % use this for pawTop
    startPosits = [40 70; 40 100; 50 70; 50 100];
    negativeEgsPerEg = 10;
    applyCircMask = true;
    circRoiPts = [55 147; 209 109; 364 136];
    negEgOverlap = .5; % allow negative examples to overlap by at most this much with the positive examples
    
    % initializations
    negImNumbers = (imNumbers(1)-1)*negativeEgsPerEg+1:imNumbers(end)*negativeEgsPerEg; % indices of negative examples
    imNumberInd = 1;
    negImNumberInd = 1;
    pixPerEg = subHgtWid(1) * subHgtWid(2);
    
    % load video and sample frame
    vid = VideoReader(file);
    wheelMask = getWheelMask(circRoiPts, [vid.Height vid.Width]);
    frame = rgb2gray(read(vid,1));
   

    % prepare figure
    close all;
    myFig = figure('units', 'normalized', 'outerposition', [0 .1 1 .9], 'keypressfcn', @keypress);
    rawAxis = subaxis(2,egsPerFrame,1:egsPerFrame, 'margin', .01', 'spacing', .01);
    rawPreview = imshow(getFeatures(frame), 'parent', rawAxis);
    circs = viscircles(gca, [0 0], 5);
    
    % for each sub frame create axes, previews, and rectangles
    for i=1:egsPerFrame
        subFrames(i).axis = subaxis(2,egsPerFrame,egsPerFrame+i);
        subFrames(i).rect = imrect(rawAxis, [startPosits(i,2) startPosits(i,1) subHgtWid(2)-1, subHgtWid(1)-1]);
        pos = getPosition(subFrames(i).rect);
        subFrames(i).img = frame(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
        subFrames(i).preview = imshow(subFrames(i).img, 'parent', subFrames(i).axis);
        addNewPositionCallback(subFrames(i).rect, @updateSubPreviews);
    end
    
    % make training set!
    while imNumberInd<=length(imNumbers)
        
        % get random frame and allow user to move rectangles around
        frame = rgb2gray(read(vid,randi(vid.numberofframes)));
        if applyCircMask
            frame = frame .* wheelMask;
        end
        set(rawPreview, 'CData', getFeatures(frame));
        updateSubPreviews();
        
        % wait for user to press enter
        stillGoing = true;
        while stillGoing; waitforbuttonpress; end
    end
    
    close all
    
    % ---------
    % FUNCTIONS
    % ---------
    
    function keypress(~,~)
        
        % save positive and create/save negative examples when ENTER is pressed
        key = double(get(myFig, 'currentcharacter'));
        
        if ~isempty(key) && isnumeric(key)
            if key==13
                stillGoing = false;
                
                % create mask of locations of positive examples
                egsMask = logical(zeros(size(frame,1), size(frame,2)));
                
                for j=1:egsPerFrame
                    pos = round(getPosition(subFrames(j).rect));
                    egsMask(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3)) = 1;
                end
                
                
                negEgPos = nan(0,2);
                % save positive examples
                for j=1:egsPerFrame
                    img = subFrames(j).img;
                    save([dataDir className '\positive\img' num2str(imNumbers(imNumberInd)) '.mat'], 'img');
                    imNumberInd = imNumberInd+1;
                    
                    % create/save negative examples for every positive example
                    for k=1:negativeEgsPerEg
                        
                        % find a frame that doesn't overlap with positive examples
                        acceptableImage = false;
                        while ~acceptableImage 
                            pos = [randi(size(frame,1)-subHgtWid(1)) randi(size(frame,2)-subHgtWid(2))]; % y,x
                            posMask = logical(zeros(size(frame,1), size(frame,2)));
                            posMask(pos(1):pos(1)+subHgtWid(1), pos(2):pos(2)+subHgtWid(2)) = 1;
                            pixelsOverlap = sum(egsMask(:)+posMask(:)>1);
                            img = frame(pos(1):pos(1)+subHgtWid(1)-1, pos(2):pos(2)+subHgtWid(2)-1,:);

                            if pixelsOverlap<(pixPerEg*negEgOverlap) && mean(img(:))>5; acceptableImage=true; end
                        end

                        % save that shit
                        save([dataDir className '\negative\img' num2str(negImNumbers(negImNumberInd)) '.mat'], 'img');
                        negImNumberInd = negImNumberInd+1;
                        negEgPos(end+1,:) = fliplr(pos+.5*fliplr(subHgtWid));
                    end
                end
                circs = viscircles(rawAxis, negEgPos, ones(1,size(negEgPos,1))*.5*mean(subHgtWid));
                waitforbuttonpress; delete(circs)
            end
        end
    end


    function updateSubPreviews(~,~)
        
        for j=1:egsPerFrame
            pos = round(getPosition(subFrames(j).rect));
            subFrames(j).img = frame(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:);
            set(subFrames(j).preview, 'CData', getFeatures(subFrames(j).img));
        end
    end
end






