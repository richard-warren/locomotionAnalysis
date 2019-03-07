close all; 


con = squeeze(controlDifs(:,1,:));
man = squeeze(modDifs(:,1,:));


figure;
plot(times, nanmean(man([flat.isBigStep],:),1)); hold on;
plot(times, nanmean(man(~[flat.isBigStep],:),1)); hold on;
plot(times, nanmean(con([flat.isBigStep],:),1));
plot(times, nanmean(con(~[flat.isBigStep],:),1));
legend('big', 'small', 'control big', 'control small');
% plot(times, nanmean(con,1));
% legend('big', 'small', 'control');





aucs = nan(2, size(con,2));
for i = 1:size(con,2)
    labels = repelem([0 1], size(con,1));
    scores = cat(1, con(:,i), man(:,i));
    
    oneBins = ~isnan(scores) & repmat([flat.isBigStep],1,2)';
    twoBins = ~isnan(scores) & repmat(~[flat.isBigStep],1,2)';
    
    if any(oneBins); [~,~,~,aucs(1,i)] = perfcurve(labels(oneBins), scores(oneBins), 1); end
    if any(twoBins); [~,~,~,aucs(2,i)] = perfcurve(labels(twoBins), scores(twoBins), 1); end
end



figure;
plot(times, aucs)
legend({'one step', 'two step'});