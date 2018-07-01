

sessions = selectSessions;
%%
isCountCorrect = nan(1,length(sessions));

for i = 1:length(sessions)

    session = sessions{i};

    vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
    locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get tracking data
    camMetaData = dlmread([getenv('OBSDATADIR') 'sessions\' session '\run.csv']); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)

    fprintf('\n%s\n', session)
    fprintf('ground truth............................ %i\n', size(camMetaData,1))
    fprintf('matlab count............................ %i\n', vid.NumberOfFrames);
    fprintf('moviepy count........................... %i\n', height(locationsTable));
%     fprintf('amount predicted by duration............ %i\n', length([0:.004:vid.Duration]));
%     fprintf('amount predicted by rounded duration.... %i\n', length([0:.004:round(vid.Duration*100)/100]));
    
    isCountCorrect(i) = size(camMetaData,1)==height(locationsTable);
end
disp('all done!')