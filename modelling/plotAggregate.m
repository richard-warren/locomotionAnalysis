%% prepare predictor

% temp
aggregate = aggregates('velocity',:);

% settings
% s.sort

% initialiaztions
img = aggregate.aggregate{1};
img = img - img(:,1);  % subtract baseline
% img = zscore(img,0,2);  % z score rows
x = linspace(aggregate.xLims(1), aggregate.xLims(2), size(img,2));
figure; imagesc(img)

%% pca
[coeff, score, ~, ~, explained] = pca(img);
figure;
pcsToShow = 5;
subplot(2,1,1); plot(cumsum(explained));
subplot(2,1,2); plot(coeff(:,1:pcsToShow), 'LineWidth', 2);

% todo: for each feature, only project portion in the meat of distribution?

%% fit gmms
maxGroups = 10;
aic = nan(1, maxGroups);
bic = nan(1, maxGroups);
gmm = cell(1, maxGroups);
pcs = 5;

for i = 1:maxGroups
    gmm{i} = fitgmdist(score(:,1:pcs), i, ...
        'RegularizationValue', .01, 'Replicates', 10, 'Start', 'randSample');
    aic(i) = gmm{i}.AIC;
    bic(i) = gmm{i}.BIC;
end

% select best number of groups
[~, nGroups] = min(bic);
groups = gmm{nGroups}.cluster(score(:,1:pcs));
fprintf('BIC groups number: %i\n', nGroups)

% plot results
figure('color', 'white', 'position', [680.00 190.00 424.00 788.00]);
subplot(2,1,1); hold on
plot(aic)
plot(bic)
legend('aic', 'bic')
title('BIC and AIC scores')

subplot(2,1,2); hold on;
groupMeans = coeff(:,1:pcs) * gmm{nGroups}.mu';  % project back onto PCs
plot(groupMeans, 'LineWidth', 2)
title('group means')




%% sort
[~, maxInd] = max(img,[],2);

% [~, sortInds] = sort(score(:,1));  % first pc projection
% [~, sortInds] = sort(aggregate.mi{1});  % mutual information
% [~, sortInds] = sort(groups);  % gmm groups
[~, sortInds] = sort(maxInd, 'descend');  % peak response time



%% plot
figure('color', 'white', 'Position', [1338.00 78.00 300.00 859.00]); hold on

imagesc(x, 1:size(img,1), img(sortInds,:))
set(gca, 'YDir', 'normal')

% add groups
groupColors = lines(nGroups);
for i = 1:nGroups
    xScat = ones(1, size(img,1)) * x(1);
    yScat = 1:size(img,1);
    yScat(groups(sortInds)~=i) = nan;
    plot(xScat, yScat, 'LineWidth', 4, 'Color', groupColors(i,:));
%     scatter(xScat, yScat, 20, groupColors(groups(sortInds),:), 'filled')
end
set(gca, 'xlim', [x(1) x(end)])


