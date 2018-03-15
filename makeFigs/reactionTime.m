

% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
neighborNum = 20;
dt = .004; % 1/fps of camera
dtInterp = .001;


% initializations
dataRaw = getKinematicData(sessions);
data = dataRaw([dataRaw.oneSwingOneStance]);
swingMaxSmps = size(data(1).modifiedLocations{1},3); % inherits max samples for swing locations from getKinematicData
times = linspace(-swingMaxSmps*dt, swingMaxSmps*dt, swingMaxSmps*2+1);
timesInterp = times(1):dtInterp:times(end);

%%

% data = dataRaw([dataRaw.oneSwingOneStance]);
% randInds = randperm(length(data), round(length(data)*.2));
% data = data(randInds);

controlVels1 = cellfun(@(x) x(1,3), {data.controlWheelVels});
controlVels2 = cellfun(@(x) x(2,3), {data.controlWheelVels});
modVels = cellfun(@(x) x(1,3), {data.modifiedWheelVels});
contactInds = reshape([data.pawObsPosInd]',4,length(data))'; contactInds = contactInds(:,3);
contactInds = contactInds + 1; % !!! this is a hack - need to fix getKinematicData s.t. no inds are zero

controlLocationsRaw1 = cellfun(@(x) x{3}(end,:,:), {data.controlLocations}, 'uniformoutput', 0);
controlLocationsRaw1 = cat(1,controlLocationsRaw1{:});
controlLocationsRaw2 = cellfun(@(x) x{3}(end,:,:), {data.controlLocations}, 'uniformoutput', 0);
controlLocationsRaw2 = cat(1,controlLocationsRaw2{:});
modLocationsRaw = cellfun(@(x) x{3}(1,:,:), {data.modifiedLocations}, 'uniformoutput', 0);
modLocationsRaw = cat(1,modLocationsRaw{:});

% subtract x and y values s.t. they start at zero at the start of every trial
controlLocationsRaw1 = controlLocationsRaw1 - repmat(controlLocationsRaw1(:,:,1), 1, 1, swingMaxSmps);
controlLocationsRaw2 = controlLocationsRaw2 - repmat(controlLocationsRaw2(:,:,1), 1, 1, swingMaxSmps);
modLocationsRaw = modLocationsRaw - repmat(modLocationsRaw(:,:,1), 1, 1, swingMaxSmps);



% initializations data containers
modLocations = nan(length(data), 2, length(times));
controlLocations = nan(length(data), 2, length(times), neighborNum);
modDifs = nan(length(data), 2, length(times));
controlDifs = nan(length(data), 2, length(times), neighborNum);


for i = 1:length(data)
    
    % get temporal inds
    trialInds = (swingMaxSmps+1 : swingMaxSmps*2) - (contactInds(i)+1);
    
    % get inds for nearest neighbor trials
    inds2 = knnsearch(controlVels2', modVels(i), 'k', neighborNum);
    
    % get control locations and dif of control locations from avg control locations
    for j = 1:length(inds2)
        
        % get single control locations
        controlLocations(i,:,trialInds,j) = controlLocationsRaw2(inds2(j),:,:);
        
        inds1 = knnsearch(controlVels1', controlVels2(inds2(j)), 'k', neighborNum);
        
        % get dif between average of remaining control locations and control trial
%         indsSub = 1:neighborNum; indsSub = indsSub(indsSub~=j);
        controlMean = nanmean(controlLocationsRaw1(inds1,:,:), 1);
        controlDifs(i,:,trialInds,j) = abs(controlLocations(i,:,trialInds,j) - controlMean);
    end
    
    % get mod locations and dif from control
    modLocations(i,:,trialInds) = modLocationsRaw(i,:,:);
    
    % get dif between average of control locations and mod locations
    controlMean = nanmean(controlLocationsRaw2(inds2,:,:), 1);
    modDifs(i,:,trialInds) = abs(modLocations(i,:,trialInds) - controlMean);    
end


% get x values
controlDifsX = squeeze(controlDifs(:,1,:,:));
controlDifsX = reshape(permute(controlDifsX, [1 3 2]), [], length(times));
modDifsX = squeeze(modDifs(:,1,:));

% interpolate
% controlDifsXInterp = interpWithNans(controlDifsX, times, timesInterp, 'pchip');
% modDifsXInterp = interpWithNans(modDifsX, times, timesInterp, 'pchip');


%% plot mean difs
xlims = [-.01 .05];

figure('color', 'white', 'menubar', 'none');

shadedErrorBar(times*1000, controlDifsX*1000, {@nanmean, @(x) nanstd(x)}, ...
    'lineprops', {'linewidth', 3, 'color', [0 0 0]}); hold on;
shadedErrorBar(times*1000, modDifsX*1000, {@nanmean, @(x) nanstd(x)}, ...
    'lineprops', {'linewidth', 3, 'color', winter(1)});
% /sqrt(size(x,1))

% pimp fig
set(gca, 'xlim', xlims*1000)
ylabel('\Deltax (mm)')
xlabel('time (ms)')
line([0 0], get(gca,'ylim'), 'color', [0 0 0])

saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes.png']);





