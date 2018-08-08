%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% get kinematic data
obsPos = -0.0087;
kinData = getKinematicData4(sessionInfo.session, obsPos);
data = kinData; save([getenv('OBSDATADIR') 'kinematicData\whiskerTrimKinematicData.mat'], 'data');

% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    preOrPost = sessionInfo.preOrPost{sessionInfoBin};
    kinData(i).preOrPost = preOrPost;
end


%% get avoidance and speed data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, false);
data = speedAvoidanceData; save([getenv('OBSDATADIR') 'kinematicData\whiskerTrimSpeedAvoidanceData.mat'], 'data');

% incorporate condition information into kinData struct
for i = 1:length(speedAvoidanceData)
    sessionInfoBin = strcmp(sessionInfo.session, speedAvoidanceData(i).session);
    preOrPost = sessionInfo.preOrPost{sessionInfoBin};
    speedAvoidanceData(i).preOrPost = preOrPost;
end

%% plot baseline kinematics

mice = unique(sessionInfo.mouse);

for j = 1:length(mice)
    
    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    preBins = strcmp(sessionInfo.preOrPost, 'pre');
    postBins = strcmp(sessionInfo.preOrPost, 'post');
    
    preSessions = sessionInfo.session(mouseBins & preBins);
    postSessions = sessionInfo.session(mouseBins & postBins);

    conditionBins = nan(1,length(kinData));
    conditionBins(ismember({kinData.session}, preSessions)) = 1;
    conditionBins(ismember({kinData.session}, postSessions)) = 2;

    conditionLabels = {'pre trimming', 'post trimming'};
    plotBaselineKinematics(kinData, conditionBins, conditionLabels, mice{j})
end


%% plot success and speed averages

% settings
circSize = 75;
lineThickness = 4;
jitterRange = .2;
lineEdges = [-.2 .2];
dvs = {'success rate', 'speed (m/s)'};

close all;
figure('color', 'white', 'menubar', 'none', 'inverthardcopy', 'off', 'position', [100 100 500 720]);

conditionNames = {'whiskers + vision', 'whiskers', 'vision', 'neither'};
wiskBins = strcmp({speedAvoidanceData.preOrPost}, 'pre');
lightOnBins = [speedAvoidanceData.isLightOn];
dvBins = {wiskBins & lightOnBins
          wiskBins & ~lightOnBins
          ~wiskBins & lightOnBins
          ~wiskBins & ~lightOnBins};

mice = unique(sessionInfo.mouse);
colors = winter(length(mice));
jitters = linspace(-jitterRange/2, jitterRange/2, length(mice));
successes = nan(length(mice), length(conditionNames));
vels = nan(length(mice), length(conditionNames));
yLims = [.4 1; .1 .8];


for i = 1:length(conditionNames)
    for mouse = 1:length(mice)

        % get session performance
        mouseBins = strcmp({speedAvoidanceData.mouse}, mice{mouse});
        
        successBins = mouseBins & dvBins{i};
        successes(mouse,i) = mean([speedAvoidanceData(successBins).isObsAvoided]);
        
        velBins = ~[speedAvoidanceData.isWheelBreak] & dvBins{i} & mouseBins;
        vels(mouse,i) = mean([speedAvoidanceData(velBins).avgVel]);

        % plot mouse avoidance and speed
        subplot(2,1,1)
        scatter(i+jitters(mouse), successes(mouse,i), circSize, colors(mouse,:), 'filled'); hold on
        subplot(2,1,2)
        scatter(i+jitters(mouse), vels(mouse,i), circSize, colors(mouse,:), 'filled'); hold on
    end

    % plot mean avoidance and speed
    subplot(2,1,1)
    line(lineEdges+i, repmat(nanmean(successes(:,i)),1,2), 'linewidth', lineThickness, 'color', 'black')
    
    subplot(2,1,2)
    line(lineEdges+i, repmat(nanmean(vels(:,i)),1,2), 'lineWidth', lineThickness, 'color', 'black')
end

% add mouse lines
for mouse = 1:length(mice)
    
    % success
    subplot(2,1,1)
    ln = line([1:length(conditionNames)] + jitters(mouse), successes(mouse,:), 'color', [.5 .5 .5]);
    uistack(ln, 'bottom')
    % vel
    subplot(2,1,2)
    ln = line([1:length(conditionNames)] + jitters(mouse), vels(mouse,:), 'color', [.5 .5 .5]);
    uistack(ln, 'bottom')
end


% pimp figs
for i = 1:2
    subplot(2,1,i)
    ylabel(dvs{i})
    set(gca, 'XTick', 1:length(conditionNames), 'XTickLabel', conditionNames, 'YLim', yLims(i,:))
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\sensoryDependencePlots.png']);
savefig([getenv('OBSDATADIR') 'figures\sensoryDependenceBarPlots.fig'])


