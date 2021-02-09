%% fit epoch glms (train in or out of reward epochs, and test generalization)


%% fit single epoch glm
session = '200622_000'; unit = 264;
[models, fitdata] = fitEpochGlm(session, unit, 'parallel', true, 'save', true);


%% fit all epoch glms

data = getUnitInfo();

tic
parfor i = 1:height(data)
    try
        fitEpochGlm(data.session{i}, data.unit(i), 'parallel', false, 'closeFig', true);
    catch exception
        fprintf('%s (%i): PROBLEM! -> %s\n', data.session{i}, data.unit(i), exception.identifier)
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)


%% collect data for plots

% settings
x = linspace(-2, 4, 200);

% get data tbl and add reponse columns
data = getUnitInfo();
m = height(data);
n = length(x);
tbl = table(nan(m,n), nan(m,n), nan(m,n), nan(m,n), nan(m,1), ...
    'VariableNames', {'resp', 'pred_in', 'pred_out', 'pred_full', 'p'});
data = cat(2, data, tbl);


% get all responses
for i = 1:height(data)
    disp(i/height(data))
    
    modelname = sprintf('E:\\lab_files\\paper2\\modelling\\glms\\epoch_glms\\%s_cell_%i_glm.mat', data.session{i}, data.unit(i));
    
    if exist(modelname, 'file')
        d = load(modelname);
        r = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'runAnalyzed', [data.session{i} '_runAnalyzed.mat']), ...
            'rewardTimes', 'omissionTimes');
        reward_all = [r.rewardTimes; r.omissionTimes];

        t        = d.fitdata.t;
        y        = d.fitdata.yRate;                    % actual firing rate
        yhatout  = d.models{'postreward', 'yhatout'};  % predicted from model trained outside rewards
        yhatin   = d.models{'postreward', 'yhatin'};   % predicted from model trained inside rewards
        yhatfull = d.models{'full', 'yhatin'};
        
        % normalize
        mn = mean(y);
        stdev = std(y);
        y        = (y        - mn) / stdev;
        yhatout  = (yhatout  - mn) / stdev;
        yhatin   = (yhatin   - mn) / stdev;
        yhatfull = (yhatfull - mn) / stdev;
%         
        % remove prediction outliers
        yhatout(yhatout<min(y))   = min(y);
        yhatout(yhatout>max(y))   = max(y);
        yhatin(yhatin<min(y))     = min(y);
        yhatin(yhatin>max(y))     = max(y);
        yhatfull(yhatfull<min(y)) = min(y);
        yhatfull(yhatfull>max(y)) = max(y);

        % compute responses for each trial
        unit_resp     = interp1(t, y, reward_all + x);         % response matrix for unit
        unit_predout  = interp1(t, yhatout, reward_all + x);   % predicted response matrix for unit
        unit_predin   = interp1(t, yhatin, reward_all + x);    % predicted response matrix for unit
        unit_predfull = interp1(t, yhatfull, reward_all + x);  % predicted response matrix for unit
        
        % add mean responses to table
        data.resp(i,:)      = nanmean(unit_resp, 1);
        data.pred_out(i,:)  = nanmean(unit_predout, 1);
        data.pred_in(i,:)   = nanmean(unit_predin, 1);
        data.pred_full(i,:) = nanmean(unit_predfull, 1);

        % response significance
        pre  = nanmean(unit_resp(:,x<=0), 2);
        post = nanmean(unit_resp(:,x>0), 2);
        [~, data.p(i)] = ttest(pre, post);
    end
end
data.delta_in   = data.resp - data.pred_in;
data.delta_out  = data.resp - data.pred_out;
data.delta_full = data.resp - data.pred_full;

%% plot reward vs predicted reward for dft models

% settings
nclusters = 5;
pcutoff = inf;
mahamax = inf;
clims = [-4 4];
% pcutoff = prctile(data.p, 80);

% restrict to modulated cells
bins = data.p < pcutoff;
data_sub = data(bins,:);
nrows = sum(bins);

close all
resp = data_sub.resp;
% resp = resp(:, x>0 & x<2) - nanmean(resp(:, x<0), 2);
% resp = resp(:, x>0 & x<2);
[clusterIds, mahaDistance, projection] = clusterResponses(resp, ...
    'plot', true, 'nclusters', nclusters, 'pcs', 5);
clusterIds(mahaDistance>mahamax) = nclusters+1;

% sort!
respAtOne = data_sub.resp(:, knnsearch(x', 1));  % activation at x=1
[sortedIds, sortInds] = sortrows([clusterIds, mahaDistance]);
sortedIds = sortedIds(:,1);

%


figure('color', 'white', 'menubar', 'none', 'position', [1.00 1.00 2560.00 1417.00])
colors = [lines(nclusters); 0 0 0];
cmap = customcolormap([0 .5 1], [1 .2 .2; 1 1 1; .2 .2 1]);

% pack data for each column into table
cols = table({data_sub.resp, data_sub.pred_full, abs(data.delta_full), data_sub.pred_out, abs(data.delta_out), data_sub.pred_in, abs(data_sub.delta_in)}', ...
             {'actual', ...
             'predicted (full)', 'abs(actual - predicted (full))', ...
             'predicted (out)', 'abs(actual - predicted (out))', ...
             'predicted (in)', 'abs(actual - predicted (in))'}', ...
             'VariableNames', {'resp', 'name'});

% without abs()
% cols = table({data_sub.resp, data_sub.pred_full, data.delta_full, data_sub.pred_out, data.delta_out, data_sub.pred_in, data_sub.delta_in}', ...
%              {'actual', ...
%              'predicted (full)', 'abs(actual - predicted (full))', ...
%              'predicted (out)', 'abs(actual - predicted (out))', ...
%              'predicted (in)', 'abs(actual - predicted (in))'}', ...
%              'VariableNames', {'resp', 'name'});

n = height(cols);

for i = 1:n
    
    resp = cols{i, 'resp'}{1};
    
    % average responses
    subplot(5, n, i); hold on;
    for j = 1:(nclusters+1)
        bins = clusterIds==j;
        mn = nanmean(resp(bins,:), 1);
        stdev = nanstd(resp(bins,:), 1);
        stdev = stdev / sqrt(sum(bins));  % (optional) convert to SEM
        patch([x fliplr(x)], [(-stdev+mn) fliplr(stdev+mn)], colors(j,:), ...
                'FaceAlpha', .25, 'EdgeColor', 'none')     % shaded error bars
        plot(x, mn, 'color', colors(j,:), 'lineWidth', 2)  % mean
    end
    if i==1; ylims = ylim; end
    plot([0 0], ylims, 'color', get(gca, 'xcolor'))     % vertical line at t=0
    title(cols{i, 'name'})
    set(gca, 'YLim', ylims)
    
    % heatmaps
    subplot(5, n, (n:n:4*n) + i); hold on
    imagesc(x, 1:nrows, resp(sortInds,:), clims)
    colormap(cmap);
    set(gca, 'ylim', [1 nrows], 'tickdir', 'out', 'ycolor', 'none')
    for j = 1:(nclusters+1)
        if any(sortedIds==j)
            ystart = find(sortedIds==j, 1, 'first')-.5;
            yend = find(sortedIds==j, 1, 'last')+.5;
            plot(x(1)*[1 1], [ystart yend], ...
                'linewidth', 3, 'color', colors(j,:))
        end
    end
    
end

saveas(gcf, 'E:\lab_files\paper2\modelling\glms\epoch_glms\summary\epoch_models.fig')


%% scatters of response magnitude before and after


% settings
clims = [-1 1];
colors = [lines(3); zeros(2,3)];  % nuclei colors
cmap = customcolormap([0 .5 1], [1 .2 .2; 1 1 1; .2 .2 1]);  % blue - white - red cmap for activation scatters

% pre - post
% prebins = x<0;
% resp_mag   = nanmean(data_sub.resp(:, ~prebins), 2) - nanmean(data_sub.resp(:, prebins), 2);
% delta_full = nanmean(data_sub.pred_full(:, ~prebins), 2) - nanmean(data_sub.pred_full(:, prebins), 2);
% delta_in   = nanmean(data_sub.pred_in(:, ~prebins), 2) - nanmean(data_sub.pred_in(:, prebins), 2);
% delta_out  = nanmean(data_sub.pred_out(:, ~prebins), 2) - nanmean(data_sub.pred_out(:, prebins), 2);

% post
prebins = x<0;
resp_mag   = nanmean(data_sub.resp(:, ~prebins), 2);
delta_full = nanmean(data_sub.pred_full(:, ~prebins), 2);
delta_in   = nanmean(data_sub.pred_in(:, ~prebins), 2);
delta_out  = nanmean(data_sub.pred_out(:, ~prebins), 2);

close all
figure('color', 'white', 'menubar', 'none', 'position', [2.00 597.00 1489.00 759.00])
lims = [min(resp_mag) max(resp_mag)];
tbl = table({delta_full, delta_in, delta_out}', ...
            {'full', 'in', 'out'}', ...
            'VariableNames', {'data', 'name'});
ccf = loadCCF();

% get scatter colors by nucleus 
[nuclei, ~, nucleiInds] = unique(data_sub.nucleus);
nucBins = ismember(data_sub.nucleus, {'fastigial', 'interpositus', 'dentate'});
scatColors = colors(nucleiInds,:);


for i = 1:height(tbl)
    pred_mag = tbl{i, 'data'}{1};  % response predicted by column model
    
    subplot(2, height(tbl), i); hold on
    plot([lims(1) lims(2)], [lims(1) lims(2)], 'color', get(gca, 'xcolor'))
    scatter(resp_mag, pred_mag, 50, scatColors, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .4)
    xlabel('response (average post)')
    ylabel('predicted response (average post)')
    set(gca, 'xlim', lims, 'ylim', lims)
        title(tbl{i, 'name'})
    
    % add legend
    if i==1
        scats = nan(1,length(nuclei));
        for j = 1:length(nuclei); scats(j) = scatter(nan, nan, [], colors(j,:), 'filled'); end
        legend(scats, nuclei, 'location', 'best')
    end
    
    % histo
    subplot(2, height(tbl), i+height(tbl)); hold on
%     temp = pred_mag; temp(temp>clims(2)) = clims(2); temp(temp<clims(1)) = clims(1);
    respColors = interp1(linspace(clims(1), clims(2), size(cmap,1))', ...
        cmap, resp_mag, 'nearest');
    plotLabels2D(ccf.labels, 'dim', 'ap', 'colors', zeros(6,3), 'patchArgs', {'facecolor', 'none'}, ...
        'apGrid', ccf.ap, 'mlGrid', ccf.ml, 'dvGrid', ccf.dv)
    scatter(data_sub.ccfMm(nucBins,1), data_sub.ccfMm(nucBins,3), ...
        [], respColors(nucBins,:), 'filled', 'MarkerFaceAlpha', .8, 'markeredgecolor', 'black')
    set(gca, 'visible', 'off')
end

saveas(gcf, 'E:\lab_files\paper2\modelling\glms\epoch_glms\summary\scatters.fig')


