% function sessionTypeCompare

% these plots test effects of dft session conditions on overall speed of mouse
% i want to know whether mice run more slowly when there are no interleaved light on trials
% also, i want to compare sessions with and without iso
% so here i plot speed as a function of obs position for sessions with and without interleaved light on trials, as well as sessions with and without markers applied after iso administration


conditions = {'all light off', 'interleaved light on', 'all light off (iso)'};
sessions = {{'180113_000', '180113_001', '180113_002', '180114_000', '180114_002', '180115_001', '180116_000', '180117_001'},...
            {'180112_000', '180112_001', '180112_002', '180114_001', '180115_000', '180115_002', '180116_001', '180117_000'},...
            {'180109_000', '180109_001', '180109_002', '180118_000', '180118_001', '180118_002'}};



% user settings
obsPrePost = [.6 .25]; % plot this much before and after the obstacle reaches the mouse
posRes = .001; % resolution of x axis, in meters
yLims = [.1 .6]; % m/s
obsPos = .382; % m, position at which obstacle is in the middle of the frame // use getFrameTimes function to determine this value

% initializations
posInterp = -obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values
totalSessions = sum(cellfun(@length, sessions));
data(totalSessions) = struct(); % stores trial data for all sessions
dataInd = 1;


% collect data
for h = 1:length(conditions)
    for i = 1:length(sessions{h})

        % load session data
        load([getenv('OBSDATADIR') 'sessions\' sessions{h}{i} '\runAnalyzed.mat'],...
                'obsPositions', 'obsTimes',...
                'obsOnTimes', 'obsOffTimes',...
                'obsLightOnTimes', 'obsLightOffTimes',...
                'wheelPositions', 'wheelTimes', 'targetFs');
        obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

        % compute velocity
        vel = getVelocity(wheelPositions, .5, targetFs);


        % iterate over all trials
        condVels = nan(length(obsOnTimes), length(posInterp));
        isLightOn = false(length(obsOnTimes), 1);
        obsOnPositions = nan(length(obsOnTimes), 1);

        for j = 1:length(obsOnTimes)

            % locate trial
            obsOnPos = obsPositions( find(obsTimes >= obsOnTimes(j), 1, 'first') );
            obsTime  = obsTimes(find( obsTimes >= obsOnTimes(j) & obsPositions >= obsPos, 1, 'first')); % time at which obstacle reaches obsPos
            obsWheelPos = wheelPositions(find(wheelTimes>=obsTime, 1, 'first')); % position of wheel at moment obstacle reaches obsPos

            % get trial positions and velocities
            trialInds = (wheelPositions > obsWheelPos-obsPrePost(1)) & (wheelPositions < obsWheelPos+obsPrePost(2));
            trialPos = wheelPositions(trialInds);
            trialPos = trialPos - obsWheelPos; % normalize s.t. 0 corresponds to the position at which the obstacle is directly over the wheel
            trialVel = vel(trialInds);

            % remove duplicate positional values
            [trialPos, uniqueInds] = unique(trialPos, 'stable');
            trialVel = trialVel(uniqueInds);

            % interpolate velocities across positional grid and store results
            trialVelInterp = interp1(trialPos, trialVel, posInterp, 'linear');
            condVels(j,:) = trialVelInterp;

            % find whether light was on
            if ~isempty(obsLightOnTimes)
                isLightOn(j) = min(abs(obsOnTimes(j) - obsLightOnTimes)) < 1; % did the light turn on near whether the obstacle turned on
            end

            % record position at which obstacle turned on
            obsOnPositions(j) = obsPos - obsOnPos;

        end
        
        % store results
        data(dataInd).lightOffVel = nanmedian(condVels(~isLightOn,:),1);
        data(dataInd).condition = conditions{h};
        data(dataInd).avgObsOnPos = nanmean(obsOnPositions);
        dataInd = dataInd + 1;
    end
end




% prepare figure
condColors = winter(2); % light on, light off
condColors = cat(1, condColors, condColors(1,:)*.4);
figure('name', 'sessionTypeCompare', 'menubar', 'none', 'units', 'pixels', 'position', [500 200 550 400], 'color', [1 1 1]);


for i = 1:length(conditions)
        
    % plot individual conditions over time
    condBins = strcmp({data.condition}, conditions{i});
    condVels = reshape([data(condBins).lightOffVel], length(posInterp), sum(condBins))';
    plot(posInterp, nanmean(condVels,1), 'color', condColors(i,:), 'linewidth', 3); hold on
        
    % plot individual iso sessions
%     if i==3
%         isoColors = copper(2);
%         isoColors = cat(1, repmat(isoColors(1,:),3,1), repmat(isoColors(2,:),3,1));
%         
%         for j = 1:size(condVels,1)
%             plot(posInterp, condVels(j,:), 'color', isoColors(j,:), 'linewidth', 1); hold on
%         end
%     end
end

% pimp fig
set(gca, 'box', 'off', 'xlim', [posInterp(1) posInterp(end)], 'ylim', yLims)
ylabel('speed (m/s)', 'fontweight', 'bold')
obsOnPos = -mean([data.avgObsOnPos]);
line([obsOnPos obsOnPos], yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2)
line([0 0], yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2)
legend(conditions)

% save fig
savefig([getenv('OBSDATADIR') 'figures\sessionTypeCompare.fig'])








