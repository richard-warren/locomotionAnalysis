

% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
neighborNum = 50;
dt = .004; % 1/fps of camera


% initializations
dataRaw = getKinematicData(sessions);
data = dataRaw([dataRaw.oneSwingOneStance]);
swingMaxSmps = size(data(1).modifiedLocations{1},3); % inherits max samples for swing locations from getKinematicData
times = linspace(-swingMaxSmps*dt, swingMaxSmps*dt, swingMaxSmps*2+1);

%%

controlVels = cellfun(@(x) x(2,3), {data.controlWheelVels});
modVels = cellfun(@(x) x(1,3), {data.modifiedWheelVels});
contactInds = reshape([data.pawObsPosInd]',4,length(data))'; contactInds = contactInds(:,3);

controlLocationsRaw = cellfun(@(x) x{3}(end,:,:), {data.controlLocations}, 'uniformoutput', 0);
controlLocationsRaw = cat(1,controlLocationsRaw{:});
modLocationsRaw = cellfun(@(x) x{3}(1,:,:), {data.modifiedLocations}, 'uniformoutput', 0);
modLocationsRaw = cat(1,modLocationsRaw{:});

% !!! subtract x and y values s.t. they start at zero at the start of every trial

% initializations data containers
modLocations = nan(length(data), 2, length(times));
controlLocations = nan(length(data), 2, length(times), neighborNum);
modDifs = nan(length(data), 2, length(times));
controlDifs = nan(length(data), 2, length(times), neighborNum);


for i = 1:length(data)
    
    % get temporal inds
    trialInds = (neighborNum+1 : neighborNum+swingMaxSmps) - (contactInds(i)+1);
    
    % get inds for nearest neighbor trials
    inds = knnsearch(controlVels', modVels(i), 'k', neighborNum+1); inds = inds(inds~=i);
    
    % get control locations and dif of control locations from avg control locations
    for j = 1:neighborNum
        
        % get single control locations
        controlLocations(i,:,j,trialInds) = controlLocationsRaw(inds(j),:,:);
        
        % get dif between average of remaining control locations and control trial
        indsSub = 1:neighborNum; indsSub = indsSub(indsSub~=j);
        controlMean = nanmean(controlLocationsRaw(indsSub,:,:), 1);
        controlDifs(i,:,trialInds,j) = abs(controlLocationsRaw(inds(j),:,:) - controlMean);
    end
    
    % get mod locations and dif from control
    modLocations(i,:,trialInds) = modLocationsRaw(i,:,:);
    
    % get dif between average of control locations and mod locations
    controlMean = nanmean(controlLocationsRaw(inds,:,:), 1);
    modDifs(i,:,trialInds) = abs(modLocationsRaw(i,:,:) - controlMean);
    
end

%%

controlDifsAvg = mean(controlDifs,4);
controlDifsX = squeeze(controlDifsAvg(:,1,:));
modDifsX = squeeze(modDifs(:,1,:));



close all; figure;
shadedErrorBar(times, controlDifsX, {@nanmean, @(x) nanstd(x)/sqrt(size(x,1))}, ...
    'lineprops', {'linewidth', 3, 'color', [.65 .65 .65]}); hold on;
shadedErrorBar(times, modDifsX, {@nanmean, @(x) nanstd(x)/sqrt(size(x,1))}, ...
    'lineprops', {'linewidth', 3, 'color', winter(1)});
set(gca, 'xlim', [times(1) times(end)])











