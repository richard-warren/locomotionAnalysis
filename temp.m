
% perform paw tracking for multiple sessions in parallel

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
views = {'bot'};
minVel = .6;
%%

for i = 1:length(sessionDirs)
    
    % get locations
    nameInd = find(sessionDirs{i}=='\',1,'last');
    try
        getLocations(sessionDirs{i}(nameInd+1:end), views, minVel, false);
    catch
        fprintf('FAILED TO ANALYZE %s\n', sessionDirs{i}(nameInd+1:end));
    end
    
end










