function locationsFixed = fixTracking(locations)

% !!! need to document // this is just post-processing for raw tracking values // interpolates short epochs of missing values and smooths! omg!

% user settings
samplesToFill = 12;
smoothSamples = 5;
interpMethod = 'spline';

% initializations
fields = fieldnames(locations);
locationsFixed = locations;

for i = 1:length(fields)
    for j = 1:size(locationsFixed.(fields{i}), 2)
        
        % fill missing values
        locationsFixed.(fields{i})(:,j) = fillShortMissing(locationsFixed.(fields{i})(:,j), samplesToFill, interpMethod);
        
        % smooth that ish
        locationsFixed.(fields{i})(:,j) = medfilt1(locationsFixed.(fields{i})(:,j), smoothSamples);
    end
end