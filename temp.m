
% perform paw tracking for multiple sessions in parallel

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');

%%

parfor i = 1:length(sessionDirs)
    
    % get locations
    nameInd = find(sessionDirs{i}=='\',1,'last');
    getLocations(sessionDirs{i}(nameInd+1:end), true);
    
end