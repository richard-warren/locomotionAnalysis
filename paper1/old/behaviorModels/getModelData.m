function data = getModelData(kinData, touchClassNames)

% takes as input kinData struct (see getKinData), and produces structure
% containing dependent variables (success, one vs two step), and predictor
% variables that can be used to build a model that predicts behavior based
% on predictor variables

% settings
touchThresh = 5;
touchesToCount = {'fore_ventral', 'hind_ventral_low'}; % only these types of touches count towards success rate



% initializations
touchInds = find(ismember(touchClassNames, touchesToCount)); % ids of touches, which appear in trialTouchesPerPaw


% copy fields that can be read directly
data = struct();
vars = {'obsPos', 'swingStartDistance', 'vel', 'angle', 'condition', 'obsHeightsVid', 'isLightOn'};
for i = 1:length(vars)
    if isfield(kinData, vars{i})
        [data(1:length(kinData)).(vars{i})] = kinData.(vars{i});
    end
end


% get remaining variables
for i = 1:length(kinData)
    modPaw = kinData(i).firstModPaw;
    
    % flip body angle w.r.t mod paw
    if modPaw==2; data(i).angle = -data(i).angle; end % angle is now wrt mod paw, but not sure in which direction...
    
    % determnine whether mouse took one big step
    data(i).oneStep = kinData(i).modStepNum(modPaw)==1; % determine
    
    % determine whether trial is successful
    data(i).totalTouchFrames = sum(any(ismember(kinData(i).trialTouchesPerPaw, touchInds), 2));
    data(i).isSuccess = data(i).totalTouchFrames<=touchThresh;
    
    % !!! need to turn 'condition' into a binary variable here
    
end
