% function prepareTrainingImages()

% temp
load('C:\Users\rick\Google Drive\columbia\obstacleData\tracking\trainingData\deepLabCut\newTrainingData\trainingData.mat', 'trainingData');
features = {'pawTL', 'pawTR', 'pawBR', 'pawBL', 'gen', 'tailBase', 'tailMid', 'tailEnd'};














% save spreadsheet for each paw, which are used by deepLabCut
for i = 1:4
    
    tableInds = (i-1)*2+2 : (i-1)*2+3;
    pawFeatures = features(:, tableInds);
    pawFeatures.Properties.VariableNames{1} = 'X';
    pawFeatures.Properties.VariableNames{2} = 'Y';
    
    writetable(pawFeatures, [writeDir 'paw' num2str(i) '.csv'], 'delimiter', ' ')
end