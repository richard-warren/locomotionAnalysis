folder = 'D:\github\chrisDlc\modifiedSpreadsheets\';

files = dir(fullfile(folder, 'originals')); files = {files(~[files.isdir]).name};

for i= 1:length(files)
    data = readtable(fullfile(folder, 'originals', files{i}));
    data = data(:,[2,1]); % flip order of two columns
    data.Properties.VariableNames = {'X', 'Y'};
    writetable(data, fullfile(folder, files{i}));
end