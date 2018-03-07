% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};


% initializations
data = getKinematicData(sessions);
dataNew = data([data.oneSwingOneStance]);





%% distance from obs

predDist = [dataNew.predictedLengths] + [dataNew.swingStartDistance];
deltaLength = cellfun(@(x) x(1,3), {dataNew.modifiedSwingLengths}) - [dataNew.predictedLengths];

% predicted distance from obs vs. swing length
close all; figure;
mice = {'run6', 'run7', 'run8'};

for i = 1:length(mice)
    subplot(length(mice),1,i)
    mouseBins = strcmp({dataNew.mouse}, mice{i});
    scatter(predDist(mouseBins), deltaLength(mouseBins))
    set(gca, 'xlim', [-.04 .03])
end

% speed distance from obs vs. swing length
% figure; scatter(cellfun(@(x) x(1,3), {dataNew.modifiedWheelVels}), deltaLength)
figure; scatter();