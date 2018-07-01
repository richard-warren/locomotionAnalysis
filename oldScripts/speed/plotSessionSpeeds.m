function plotSessionSpeeds(session)


% settings
obsPrePost = [.6 .25]; % plot this much before and after the obstacle reaches the mouse
posRes = .001; % resolution of x axis, in meters
yLims = [.1 .6]; % m/s
obsPos = .382; % m, position at which obstacle is in the middle of the frame // use getFrameTimes function to determine this value
traces = 20;
wiskPos = .330; % temp, this is an approximation
avgVelPositions = [.336 .436]; % compute trial avg velocity only between these positions


% initializations
posInterp = -obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values
avgVelPositions = avgVelPositions - obsPos; % compute trial avg velocity only between these positions

% get session info
[sessionVels, isLightOn, obsOnPositions] = getSessionSpeedInfo(session, posInterp, obsPos);
sessionVels = sessionVels(~isLightOn, :);


% plot
figure('name', [session 'TrialSpeeds'], 'color', 'white', 'position', [200 200 1000 400]);
cmap = winter(size(sessionVels,1));
traceInds = floor(linspace(1,size(sessionVels,1),traces));

subplot(1,3,2:3)
for i = traceInds
    plot(posInterp, sessionVels(i,:), 'linewidth', 1, 'color', cmap(i,:)); hold on    
end
plot(posInterp, nanmean(sessionVels,1), 'linewidth', 3, 'color', mean(cmap,1))

% pimp fig
set(gca, 'box', 'off', 'xlim', [posInterp(1) posInterp(end)], 'ylim', yLims)
xlabel('position (m)', 'fontweight', 'bold');
obsOnPos = -mean(obsOnPositions);   
line([obsOnPos obsOnPos], yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2)
line([0 0], yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2)
line(repmat(wiskPos-obsPos,1,2), yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2)


subplot(1,3,1)
velInds = posInterp>=avgVelPositions(1) & posInterp<=avgVelPositions(2);
sessionMeans = nanmean(sessionVels(:, velInds), 2);
scatter(zeros(1,size(sessionVels,1)), sessionMeans, 50, cmap, 'filled', 'jitter', 'on');
set(gca, 'ylim', yLims, 'xcolor', 'none')
ylabel('speed (m/s)', 'fontweight', 'bold')
