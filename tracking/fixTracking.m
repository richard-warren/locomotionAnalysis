function locationsFixed = fixTracking(locations)

% !!! need to document // this is just post-processing for raw tracking values // interpolates short epochs of missing values and smooths! omg!

% user settings
samplesToFill = 5;
smoothSamples = 3;
interpMethod = 'linear';

% initializations
fields = fieldnames(locations);
locationsFixed = locations;

for i = 1:length(fields)
    for j = 1:size(locationsFixed.(fields{i}), 2)
        
        % fill missing values
        locationsFixed.(fields{i})(:,j) = fillShortMissing(locationsFixed.(fields{i})(:,j), samplesToFill, interpMethod);

        % smooth that ish
        locationsFixed.(fields{i})(:,j) = medfilt1(locationsFixed.(fields{i})(:,j), smoothSamples, 'omitnan');

%         nanInds = isnan(locationsFixed.(fields{i})(:,j));
%         locationsFixed.(fields{i})(:,j) = smooth(locationsFixed.(fields{i})(:,j), smoothSamples);
%         locationsFixed.(fields{i})(nanInds,j) = nan;
    end
end