function findCcfSections(mouse)
% to perform histo registration to allen brain common coordinate framework
% (ccf), it is necessary to find the sections in ccf that correspond to the
% final manually traced sections of the nuclei in the histology // launch
% this gui to flip through ccf sections to find the one that matches the
% final traced sections of the nuclei for a given mouse

% todo: flip left right imgs // store everything in properties // mask
% nuclei

% inits
fig = figure('name', mouse, 'units', 'pixels', 'position', [61.00 455.00 1824.00 413.00],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
ccf = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));
ccfSection = 450;

% get final sections whether left and right nuclei are labelled
leftInds = find(contains(data.nuclei, 'Left'));
rightInds = find(contains(data.nuclei, 'Right'));
leftLabels = ismember(data.labels, uint8(leftInds));
rightLabels = ismember(data.labels, uint8(rightInds));
finalLeftInd = find(any(leftLabels, [2 3]), 1, 'last');
finalRightInd = find(any(rightLabels, [2 3]), 1, 'last');

% final section
isLeft = true;
updateSection();

% ccf
updateCcf(0);



% keypress controls
function changeFrames(~,~)
    
    key = double(get(fig, 'currentcharacter'));
    
    if ~isempty(key) && isnumeric(key)
        
        % left: move frame backward
        if key==28
            updateCcf(-1);
        
        % right: move frame forward
        elseif key==29                  
            updateCcf(1);
        
        % f: flip section left/right
        elseif key==102
            updateSection();
        
        % escape: close window
        elseif key==27
            close(fig)
        end
    end
end

function updateCcf(deltaSection)
    subplot(1,2,2)
    ccfSection = ccfSection + deltaSection;
    ccfSection = max(min(ccfSection, size(ccf.imgs,1)),1);
    
    mask = squeeze(ccf.labels(ccfSection,:,:))>0;
    
    img = squeeze(ccf.imgs(ccfSection,:,:));
    img(mask) = 255;
    image(img); colormap gray; axis equal
    set(gca, 'color', 'black', 'YDir', 'reverse', 'xcolor', 'none', 'ycolor', 'none')
    title(sprintf('ccf (%i/%i)', ccfSection, size(ccf.imgs,1)), 'color', 'white')
end

function updateSection()
    subplot(1,2,1);
    isLeft = ~isLeft;
    if isLeft
        ind = finalLeftInd;
        string = 'left';
    else
        ind = finalRightInd;
        string = 'right';
    end
    
    mask = squeeze(data.labels(ind,:,:))>0;
    img = squeeze(data.imgs(ind,:,:));
    img(mask) = 255;
    image(img); colormap gray; axis equal
    set(gca, 'color', 'black', 'YDir', 'reverse', 'xcolor', 'none', 'ycolor', 'none')
    title(sprintf('final %s label (%i/%i)', string, ind, size(data.imgs,1)), 'color', 'white')
end

end