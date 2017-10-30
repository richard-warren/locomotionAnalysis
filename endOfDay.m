
% preliminary spike analysis
spikeAnalysis('C:\Users\Rick\Google Drive\columbia\obstacleData\sessions')


% generate plots
obsAvoidance2('run3', {'obsTest', 'obsTestBr'})
obsAvoidance2('run4', 'obsTest')
obsAvoidance2('run5', {'obsTest', 'obsTestBr'})

% make videos
makeVid('171029_001');
