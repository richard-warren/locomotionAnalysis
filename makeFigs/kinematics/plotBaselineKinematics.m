function plotBaselineKinematics(data, conditionBins, conditionLabels, figTitle)

% temp
% conditionBins = ones(1,length(data));
% sessions = unique({data.session});
% for i = 2:length(sessions); conditionBins(strcmp({data.session}, sessions{i})) = i; end
% conditionLabels = sessions;
% conditionLabels = {'test'};

% settings
xLims = [0, .12];
yLims = [-.025 .025];
zLims = [-.008 .006];
leftDarkening = .5;

% initializations
conditions = unique(conditionBins(~isnan(conditionBins)));
colors = hsv(length(conditionLabels));

noObsLocations = cell(1, length(conditionLabels)); % one cell entry per condition
xOffsets = cell(1,length(conditionLabels));
for i = conditions
    conditionLocations = cellfun(@(x) cat(4,x{:}), {data(conditionBins==i).noObsLocationsInterp}, 'UniformOutput', 0);
    conditionLocations = cat(1, conditionLocations{:});
    
    % get avg x offset between fore and hind paws
    meanHindX = mean(mean(conditionLocations(:,1,1,[1, 4])));
    meanForeX = mean(mean(conditionLocations(:,1,1,[2, 3])));
    xOffsets{i} = meanForeX - meanHindX;
    
    % normalize s.t. hinds start on avg at 0 and fores start on avg 0 + xOffsets{i}
    for j = 1:4
        xStarts = squeeze(conditionLocations(:,1,1,j));
        conditionLocations(:,1,:,j) = squeeze(conditionLocations(:,1,:,j)) - repmat(xStarts,1,size(conditionLocations,3));
    end
    
    noObsLocations{i} = cat(1, conditionLocations);
end


% prepare figure
figure('Name', figTitle, 'color', 'white', 'menubar', 'none', 'position', [300 100 1200 800]);

% plot xy
for i = conditions
    
    pawXOffsets = [0 xOffsets{i} xOffsets{i} 0]; % offset forepaws only
    
    % get condition averages
    stepAvgs = squeeze(mean(noObsLocations{i},1));
    
    for j = 1:4
        % plot xy
        subplot(length(noObsLocations),2,1:2:length(noObsLocations)*2)
        if j==1 ||j==2; darkening = leftDarkening; else; darkening=1; end
        plot(stepAvgs(2,:,j), stepAvgs(1,:,j)+pawXOffsets(j), ...
            'linewidth', 3, 'color', colors(i,:)*darkening); hold on;
    end
    
    for j = 1:4
        if j==1 ||j==2; darkening = leftDarkening; else; darkening=1; end
        subplot(length(noObsLocations),2,2+(i-1)*2)
        plot(stepAvgs(1,:,j)+pawXOffsets(j), stepAvgs(3,:,j), ...
            'linewidth', 3, 'color', colors(i,:)*darkening); hold on;
    end
end

% pimp fig
subplot(length(noObsLocations),2,1:2:length(noObsLocations)*2)
daspect([1 1 1]);
set(gca, 'XLim', yLims, 'YLim', xLims, 'Box', 'off')

for i = 1:length(noObsLocations)
    subplot(length(noObsLocations),2,2+(i-1)*2)
    daspect([1 1 1]);
    set(gca, 'XLim', xLims, 'YLim', zLims, 'box', 'off')
    title(conditionLabels{i})
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\baselineKinematics.png']);
savefig([getenv('OBSDATADIR') 'figures\baselineKinematics.fig'])


