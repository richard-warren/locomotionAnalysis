function [a] = loadPickle(filename)
    
    % settings
    pythonPath = 'C:\\Users\\rick\\Anaconda3\\envs\\deepLabCut\\python';
    
    if ~exist(filename,'file')
        error('%s is not a file',filename);
    end
    
    outname = [tempname() '.mat'];
    pyscript = ['import pickle as pickle;import sys;import scipy.io;file=open(""' filename '"",""rb"");dat=pickle.load(file);file.close();scipy.io.savemat(""' outname '"",dat)'];
    system([pythonPath ' -c "' pyscript '"']);
    a = load(outname);
end