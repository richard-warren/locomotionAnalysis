

close all; figure;
plot(locations.locationsRaw(:,1,1), 'linewidth', 3);
hold on; plot(locations.locationsCorrected(:,1,1))

correctedInds = find(~isnan(locations.locationCorrections(:,1,4)));

hold on; scatter(correctedInds, ones(1,length(correctedInds))*140)