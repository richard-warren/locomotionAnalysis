


censoredPeriod = .001;
refractoryPeriod = .003;
time = 1000;
violations = 20;
totalSpikes = 10000;

rhs = 2 * (refractoryPeriod-censoredPeriod) * totalSpikes^2 / violations / (time); % what the fuck is this math about??? is this really correct???
disp(.5 - .5*sqrt((rhs-4)/rhs));