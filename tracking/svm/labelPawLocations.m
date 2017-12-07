function labelPawLocations(vidFile, frameInds, totalEgs)

% !!! need to document


% settings
objectNum = 4;
colors = [1 76 153; 1 153 1; 255 255 51; 255 1 255] / 255;
figSize = 2;

% initializations
vid = VideoReader(vidFile);
frame = read(vid,1);
egInd = 1;
stillGoing = true;
fileDir = vidFile(1:find(vidFile=='\',1,'last'));
fileName = vidFile(find(vidFile=='\',1,'last')+1 : find(vidFile=='.',1,'last')-1);
keyboard
frameInd = [];


% prepare figure
close all;
fig = figure('units', 'pixels', 'outerposition', [300 300 vid.Width*figSize vid.Height*figSize],...
             'color', [0 0 0], 'menubar', 'none', 'keypressfcn', @keypress);
imPreview = image(getFeatures(frame));
hold on; scatter([1, vid.Width, vid.Width, 1], [1, 1, vid.Height, vid.Height], 200, colors, 'filled');
set(gca, 'units', 'normalized', 'position', [0 0 1 1])
getNewFrame;

% initialize draggable objects
impoints = cell(1,objectNum);
locations = nan(2, totalEgs, objectNum); % ([x,y], egNum, obNum)
locationFrameInds = nan(1, totalEgs);

for i = 1:objectNum    
    impoints{i} = impoint(gca, [10 10]*i);
    setColor(impoints{i}, colors(i,:));
end


while stillGoing
    waitforbuttonpress
end

% save and close up shop
saveData();
close(fig);



% ---------
% FUNCTIONS
% ---------

function keypress(~,~)
        
    % save positive and create/save negative examples when ENTER is pressed
    key = double(get(fig, 'currentcharacter'));

    if ~isempty(key) && isnumeric(key)
        switch key
        
            % commit locations for frame
            case 13 % enter
                
                % sore locations
                for j = 1:objectNum    
                    xy = getPosition(impoints{j});
                    locations(1:2, egInd, j) = xy;
                end
                
                % store frame ind
                locationFrameInds(egInd) = frameInd;
                
                % get next frame
                getNewFrame()
                
                % update counter
                egInd = egInd + 1;
                if egInd > totalEgs
                    stillGoing = false;
                else
                    disp(egInd)
                end
            

            % get new frame
            case 110 % 'n'
                getNewFrame();
            
            % save current progress
            case 115 % 's'
                saveData();
        end
    end
end


function getNewFrame
    frameInd = frameInds(randi(length(frameInds)));
    frame = read(vid, frameInd);
    set(imPreview, 'CData', frame);
end

function saveData
    save([fileDir fileName 'handLabeledLocations'], 'locations', 'locationFrameInds');
end

end