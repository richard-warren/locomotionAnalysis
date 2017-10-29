function deleteFiles(filesToDelete)

% goes through folders in dataDir and deletes all files contained within filesToDelete
%
% input      filesToDelete:  cell array of file names to be deleted from folders within dataDir



% settings
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';


% make sure input is cell array
if ischar(filesToDelete)
    filesToDelete = {filesToDelete};
end



% get data folders
dataFolders = dir(dataDir);
dataFolders = dataFolders(3:end); % remove current and parent directory entries
dataFolders = dataFolders([dataFolders.isdir]); % keep only folders



% iternate through all folders
for i = 1:length(dataFolders)
    
    % get files in folder
    filesInDir = dir([dataDir dataFolders(i).name]);
    fileNames = {filesInDir.name};
    
    % delete all files in fileToDelete
    for j = 1:length(filesToDelete)
        
        if any(strcmp(fileNames, filesToDelete{j}))
            
            fprintf('%s: deleting %s\n', dataFolders(i).name, filesToDelete{j});
            delete([dataDir dataFolders(i).name '\' filesToDelete{j}]); 
        end 
    end
end


