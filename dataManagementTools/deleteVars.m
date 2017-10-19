function deleteVars(varsToDelete)

% deletes variable(s) contained in all run.mat files

dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';

dataFolders = dir(dataDir);
dataFolders = dataFolders(3:end); % remove current and parent directory entries
dataFolders = dataFolders([dataFolders.isdir]); % keep only folders


for i = 1:length(dataFolders)
    
    % load data
    file = [dataDir dataFolders(i).name '\run.mat'];
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


