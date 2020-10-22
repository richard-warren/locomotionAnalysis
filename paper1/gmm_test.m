% settings
componentOffset = 3;
samples = 1000;

% generate data
oneComponent = normrnd(0,1,1,samples);
twoComponent = [normrnd(0,1,1,round(samples/2)) normrnd(componentOffset,1,1,round(samples/2))];

% fit 1 and 2 component models
args = {'Replicates', 10, 'options', statset('MaxIter', 1000)};
oneComponent1 = fitgmdist(oneComponent', 1, args{:});
oneComponent2 = fitgmdist(oneComponent', 2, args{:});
twoComponent1 = fitgmdist(twoComponent', 1, args{:});
twoComponent2 = fitgmdist(twoComponent', 2, args{:});

% compute likelihood ratio test
oneComponentRatio = -2 * (twoComponent2.NegativeLogLikelihood - twoComponent1.NegativeLogLikelihood);
twoComponentRatio = -2 * (oneComponent2.NegativeLogLikelihood - oneComponent1.NegativeLogLikelihood);
fprintf('------\n')
fprintf('two-component likelihood ratio: %.2f\n', oneComponentRatio)
fprintf('one-component likelihood ratio: %.2f\n', twoComponentRatio)

% plot
figure('color', 'white', 'position', [363.00 706.00 560.00 420.00]); hold on
histo = histogram(twoComponent, 50);
histogram(oneComponent, 'BinEdges', histo.BinEdges)
xLims = xlim;
legend('two component', 'one component')