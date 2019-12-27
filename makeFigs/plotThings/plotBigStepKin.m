function plotBigStepKin(kin, kinCtl, obsHgts, rowInds, isBigStep, varargin)

% plots kinematics of first modified paw, with one line showing
% trajectories for big steps, and another for shortened steps // relative
% thickness of lines represents portion of trials on which one vs the other
% strategy is chosen // also shows histograms of showing prob of landing in
% dft spots relative to obs // kin is trialNum X 2(xz) X N, as is ctlKin // rowInds
% tells which row each trial belongs to, and isBigStep is whether each
% trial is a big step or not // obsHgts is height of obstalce for each
% trial

% to do: add individual trials option // *add times of whisker contact


% global settings
s.obsRadius = 3.175/1000/2; % (m)
s.colors = [.25 1 1; .25 1 .25]; % colors for little vs big steps
s.obsColor = [.2 .2 .2];
s.controlColor = [.4 .4 .4];
s.xLims = [-.1 .04];
s.yLims = [0 .015];
s.lineWid = 2;
s.addHistos = true;
s.histoFillAlpha = .4;
s.histoHgt = .5; % expressed as fraction of kinematic plots
s.contactInds = []; % inds for each trial at which wisk contacts obs
s.showSmpNum = false; % if true, shows number of big and small steps per bin
s.ctlDistances = []; % if provided, uses this as control distro instead of the final x position of the control kinematics for each step
s.verticalBuffer = .05;  % space to added between rows of plot, expressed as fraction of total height



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings contained in opts
if s.addHistos; s.yLims(1) = -s.yLims(2); end
rowInds=rowInds(:); isBigStep=isBigStep(:); % make sure everybody is a column vector
numRows = max(rowInds);
xGrid = linspace(s.xLims(1), s.xLims(2), 500);
kinCtl(:,1,:) = kinCtl(:,1,:) - kinCtl(:,1,1) + kin(:,1,1); % shift ctrl locations so they start at the same x positions as mod locations




for h = 1:numRows

    % get subplot bins
    axes('position', [0 1-h*(1/numRows)+s.verticalBuffer 1 1/numRows-s.verticalBuffer])
    hold on

    % get subplot bins for different conditions
    ctlBins = rowInds==h;
    oneStepBins = rowInds==h & isBigStep;
    twoStepBins = rowInds==h & ~isBigStep;
    oneTwoRatio = sum(oneStepBins) / (sum(oneStepBins) + sum(twoStepBins)); % ratio of trials in which swing foot takes one large step to those in which an additional step is taken

    % add line to bottom
    line(s.xLims, [0 0], 'color', [1 1 1]*.15, 'linewidth', 1)
    
    % plot control steps
    x = squeeze(nanmean(kinCtl(ctlBins,1,:),1));
    z = squeeze(nanmean(kinCtl(ctlBins,2,:),1));
    plot(x, z, 'color', s.controlColor, 'linewidth', s.lineWid);
    
    % plot modified steps
    allBins = {twoStepBins, oneStepBins};
    ratios = [1-oneTwoRatio, oneTwoRatio];

    for i = 1:2
        if any(allBins{i})
            % plot traces
            x = squeeze(nanmean(kin(allBins{i},1,:),1));
            z = squeeze(nanmean(kin(allBins{i},2,:),1));
            plot(x, z, 'color', s.colors(i,:), 'linewidth', s.lineWid*ratios(i)*2); hold on;
            
            if ~isempty(s.contactInds) % !!! finding mean contact ind is not actually fair, bc excludes trials in which contact ind does not occur during swing // this is a little misleading...
                ind = round(nanmean(s.contactInds(allBins{i}))); % mean contact ind
                scatter(nanmean(x(ind)), nanmean(z(ind)), 100*ratios(i), s.colors(i,:), 'filled');
            end
        end
    end

    % add obstacle
    avgObsHgt = nanmean(obsHgts(rowInds==h)); % get avg obstacle height for row
    rectangle('position', [0-s.obsRadius, avgObsHgt-2*s.obsRadius, 2*s.obsRadius, 2*s.obsRadius], ...
        'curvature', [1 1], 'facecolor', s.obsColor, 'edgecolor', 'none');
    
    % add sample size
    if s.showSmpNum
        text(s.xLims(1), s.yLims(2), ...
            sprintf('%i (%i,%i)', sum([oneStepBins; twoStepBins]), sum(oneStepBins), sum(twoStepBins)), 'VerticalAlignment', 'top')
    end
    
    % add god damn histograms
    if s.addHistos
        
        % control
        if isempty(s.ctlDistances)
            kd = ksdensity(kinCtl(ctlBins,1,end), xGrid);
        else
            kd = ksdensity(s.ctlDistances(ctlBins), xGrid);
        end
        kd = -kd * (s.yLims(2)/max(kd)) * s.histoHgt;
        fill([xGrid xGrid(end) xGrid(1)], [kd 0 kd(1)], s.controlColor, 'FaceAlpha', s.histoFillAlpha)
        plot(xGrid, kd, 'Color', s.controlColor, 'LineWidth', 2)
        
        % mod
        for i = 1:2
            if any(allBins{i})
                kd = ksdensity(kin(allBins{i},1,end), xGrid);
                kd = -kd * (s.yLims(2)/max(kd)) * s.histoHgt * ratios(i);
                fill([xGrid xGrid(end) xGrid(1)], [kd 0 kd(1)], s.colors(i,:), 'FaceAlpha', s.histoFillAlpha)
                plot(xGrid, kd, 'Color', s.colors(i,:), 'LineWidth', 2)
            end
        end
    end

    % set appearance
    set(gca, 'xlim', s.xLims, 'ylim', s.yLims, 'DataAspectRatio', [1 1 1], 'Visible', 'off');
end






