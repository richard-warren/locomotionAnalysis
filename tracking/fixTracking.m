function locationsFixed = fixTracking(locations)

% !!! need to document // this is just post-processing for raw tracking values // interpolates short epochs of missing values and smooths! omg!

% setting
interpMethod = 'spline';

% initializations
trials = unique(locations.trialIdentities(~isnan(locations.trialIdentities)));
locationsFixed = locations;


for i = trials'
    
    trialInds = find(locationsFixed.trialIdentities==i);
    
    for dimension = 1:2
        for paw = 1:4
            locationsFixed.locationsRaw(trialInds, dimension, paw) = ...
                fillmissing(locationsFixed.locationsRaw(trialInds, dimension, paw), interpMethod, 'endvalues', 'none');
        end
    end
end




% % user settings
% samplesToFill = 10;
% smoothSamples = 3;
% interpMethod = 'spline';
% 
% % initializations
% locationsFixed = locations;
% 
% for i = 1:2
%     for j = 1:4
%         
%         % fill missing values
%         locationsFixed(:,i,j) = fillShortMissing(locationsFixed(:,i,j), samplesToFill, interpMethod);
% 
%         % smooth that ish
%         locationsFixed(:,i,j) = medfilt1(locationsFixed(:,i,j), smoothSamples, 'omitnan');
% 
% %         nanInds = isnan(locationsFixed.(fields{i})(:,j));
% %         locationsFixed.(fields{i})(:,j) = smooth(locationsFixed.(fields{i})(:,j), smoothSamples);
% %         locationsFixed.(fields{i})(nanInds,j) = nan;
%     end
% end