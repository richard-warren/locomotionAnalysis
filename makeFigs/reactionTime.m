

% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
neighborNum = 50;

% initializations
dataRaw = getKinematicData(sessions);
data = dataRaw([dataRaw.oneSwingOneStance]);


for i = 1:length(data)
    
    inds = knnsearch([data.vel]', data(1).vel, 'k', neighborNum+1);
    inds = inds(inds~=i);
    
end



% swingInds = reshape([data.pawObsPosInd]',4,length(data))';

