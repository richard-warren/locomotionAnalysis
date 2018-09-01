%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'lesionNotes');
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & strcmp(sessionInfo.brainRegion, 'mtc'),:);



%% compute kinematic data

loadPreviousData = false;

if loadPreviousData
    load([getenv('OBSDATADIR') 'matlabData\lesionKinematicData.mat'], 'data');
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, data);
else
    kinData = getKinematicData4(sessionInfo.session, sessionInfo, []);
end
data = kinData; save([getenv('OBSDATADIR') 'matlabData\lesionKinematicData.mat'], 'data', '-v7.3', '-nocompression'); clear data;



%% load kinematic data

load([getenv('OBSDATADIR') 'matlabData\lesionKinematicData.mat'], 'data');
kinData = data; clear data;
disp('lesion kinematic data loaded!')


%% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    sideOfBrain = sessionInfo.side{sessionInfoBin};
    condition = sessionInfo.preOrPost{sessionInfoBin};
    
    % get which number session post lesion
    mouse = sessionInfo.mouse{sessionInfoBin};
    mousePostSessions = sessionInfo.session(strcmp(sessionInfo.mouse, mouse) & strcmp(sessionInfo.preOrPost, 'post'));
    postSessionNum = find(strcmp(mousePostSessions, sessionInfo.session{sessionInfoBin}));
    
    if strcmp(sideOfBrain, 'left'); contraLimb = 3;
    elseif strcmp(sideOfBrain, 'right'); contraLimb = 2;
    else; contraLimb = nan; end
    
    kinData(i).contraPawFirst = kinData(i).firstPawOver==contraLimb;
    kinData(i).ipsiPawFirst = kinData(i).firstPawOver~=contraLimb && ~isnan(contraLimb);
    kinData(i).brainRegion = brainRegion;
    kinData(i).sideOfBrain = sideOfBrain;
    kinData(i).condition = condition;
    kinData(i).postSessionNum = postSessionNum;
end

%% get speed and avoidance data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, false);

% incorporate condition information into kinData struct
for i = 1:length(speedAvoidanceData)
    sessionInfoBin = strcmp(sessionInfo.session, speedAvoidanceData(i).session);
    
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    sideOfBrain = sessionInfo.side{sessionInfoBin};
    condition = sessionInfo.preOrPost{sessionInfoBin};
    
    speedAvoidanceData(i).brainRegion = brainRegion;
    speedAvoidanceData(i).sideOfBrain = sideOfBrain;
    speedAvoidanceData(i).condition = condition;
end

data = speedAvoidanceData; save([getenv('OBSDATADIR') 'matlabData\lesionSpeedAvoidanceData.mat'], 'data'); clear data;


%% load speed and avoidance data

load([getenv('OBSDATADIR') 'matlabData\lesionSpeedAvoidanceData.mat'], 'data');
speedAvoidanceData = data; clear data;
disp('lesion speed avoidance data loaded!')

%% plot baseline kinematics

mice = unique(sessionInfo.mouse);

for j = 1:length(mice)
    
    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    brainRegion = unique(sessionInfo.brainRegion(mouseBins));
    preBins = strcmp(sessionInfo.preOrPost, 'pre');
    postBins = strcmp(sessionInfo.preOrPost, 'post');
    sides = unique(sessionInfo.side(mouseBins));
    
    for k = 1:length(sides)
        sideBins = strcmp(sessionInfo.side, sides{k});
        preSessions = sessionInfo.session(mouseBins & sideBins & preBins);
        postSessions = sessionInfo.session(mouseBins & sideBins & postBins);
        
        conditionBins = nan(1,length(kinData));
        conditionBins(ismember({kinData.session}, preSessions)) = 1;
        conditionBins(ismember({kinData.session}, postSessions)) = 2;
        
        conditionLabels = {['pre - ' sides{k} ' ' brainRegion{1}], ['post - ' sides{k} ' ' brainRegion{1}]};
        plotBaselineKinematics(kinData, conditionBins, conditionLabels, mice{j})
    end
end


%% plot dv averages

manipulationBarPlots(speedAvoidanceData, {'pre', 'post'}, 'lesions');
saveas(gcf, [getenv('OBSDATADIR') 'figures\lesions\lesionBarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\lesions\leionsBarPlots.fig'])


%% plot dvs across sessions

manipulationAcrossSessions(speedAvoidanceData, {'pre', 'post', 'postNoWisk'}, 'lesions');
saveas(gcf, [getenv('OBSDATADIR') 'figures\lesions\lesionAcrossSessions.png']);
savefig([getenv('OBSDATADIR') 'figures\lesions\lesionAcrossSessions.fig'])


%% paw height by obs height for mtc

binNames = {'ipsi', 'contra', 'pre'};

bins = zeros(1,length(kinData));
bins(strcmp({kinData.condition}, 'post') & strcmp({kinData.brainRegion}, 'mtc') & [kinData.ipsiPawFirst]) = 1;
bins(strcmp({kinData.condition}, 'post') & strcmp({kinData.brainRegion}, 'mtc') & [kinData.contraPawFirst]) = 2;
bins(strcmp({kinData.condition}, 'pre')) = 3;

scatterObsVsPawHeights(kinData, bins, binNames);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/lesions/heightShapingScatter.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/lesions/heightShapingScatter.fig'))


%% plot kinematics

% get trial bins
contraFirstBins = [kinData.contraPawFirst];
ipsiFirstBins = [kinData.ipsiPawFirst];
preBins = strcmp({kinData.condition}, 'pre');
postBins = strcmp({kinData.condition}, 'post');
lightOnBins = [kinData.isLightOn];
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');

% mtc
bins = zeros(1,length(kinData));
bins(mtcBins & ipsiFirstBins & postBins) = 1; % ipsi
bins(mtcBins & contraFirstBins & postBins) = 2; % contra
bins(mtcBins & preBins) = 3; % pre
binLabels = {'ipsi', 'contra', 'pre'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures/lesions/lesionKinematics.png'));
savefig(fullfile(getenv('OBSDATADIR'), 'figures/lesions/lesionKinematics.fig'))


%% edit example videos for talks

dir = 'C:\Users\rick\Desktop\talkVids\';
startTrial = 150;
trialsToShow = 10;
iti = 5;
trials = startTrial:iti:startTrial+trialsToShow*5;
makeVidWisk(fullfile(dir, 'preLesion'), '180808_000', [-.05 .1], .1, trials');
makeVidWisk(fullfile(dir, 'postLesion'), '180811_000', [-.05 .1], .1, trials');






