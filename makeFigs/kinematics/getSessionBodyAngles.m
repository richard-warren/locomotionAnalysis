function angles = getSessionBodyAngles(locationsTable, nosePos)

% computes the angle between the nose and the base of the tail for every
% frame // positive angles are leaning towards left side of mouse body, and negative
% angles are leaning towards right side


% load kinematic data
% load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'nosePos')
% try
%     locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
%     tailXY = [locationsTable.tailBase_bot, locationsTable.tailBase_bot_1];
% catch
%     disp('hacky fix') % !!! for some reason readtable sometimes fails, and reading first with xls read solves the problem... wtf
%     locationsTable = xlsread([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
%     locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
%     tailXY = [locationsTable.tailBase_bot, locationsTable.tailBase_bot_1];
% end

tailXY = [locationsTable.tailBase_bot, locationsTable.tailBase_bot_1];
tailXY(:,1) = -(tailXY(:,1) - nosePos(1)); % set X to number of pixels behind nose
tailXY(:,2) = -(tailXY(:,2) - nosePos(2));  % set y to number of pixels above nose

angles = atan2(tailXY(:,2), tailXY(:,1));
angles = rad2deg(angles);