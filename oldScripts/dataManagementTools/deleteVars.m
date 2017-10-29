function deleteVars(varsToDelete, matFile)

% goes through all folders in dataDir, and deletes the variables listed in cell array varsToDelete contained within matFile
%
% input      varsToDelete:  cell array of variable names to be deleted from mat file
%            matFile:       name of a .mat file from which variables will be deleted (e.g. 'run.mat')

dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';

dataFolders = dir(dataDir);
dataFolders = dataFolders(3:end); % remove current and parent directory entries
dataFolders = dataFolders([dataFolders.isdir]); % keep only folders


for i = 1:length(dataFolders)
    
    % load data
    file = [dataDir dataFolders(i).name '\' matFile];
    data = load(file);
    deleted = false;
    
    % delete variables
    for j = 1:length(varsToDelete)
        
        if any(strcmp(fieldnames(data), varsToDelete{j}))
            
            fprintf(['deleting ''' varsToDelete{j} ''' in %s\n'], dataFolders(i).name)
            data = eval(['rmfield(data, ''' varsToDelete{j} ''')']);
            deleted = true;
            
        end
    end
    
    % save file if anything was deleted
    if deleted
            save(file, '-struct', 'data');
    else
        fprintf('variables not found in %s\n', dataFolders(i).name)
    end
    
end


