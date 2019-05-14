%% COMPUTE PREDICTORS

% global settings
referenceModPaw = 2;
velSmps = 10; % how many samples to compute velocity over
frameRate = 250;

% load experiment data
fprintf('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

%% flatten data and compute predictors
flat = flattenData(data, {'mouse', 'session', 'trial', ...
    'firstModPaw', 'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'isBigStep'});
flat = struct2table(flat);

% flip predictors in flat so everything is relative to firstModPaw
flipBins = [flat.firstModPaw]==referenceModPaw;
flat.angleAtWiskContact(flipBins) = -flat.angleAtWiskContact(flipBins);


%% SET UP PREDICTORS


sessions = unique(flat.session);

for i = 1:length(sessions)
    
    fprintf('%s: computing extra predictor variables\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'kinData')
    sessionBins = strcmp(flat.session, sessions{i});
    if sum(sessionBins)~=length(kinData); disp('WTF! mismatch in number of trials! fucking shit!!!'); end
    
    % !!! make sure 0s don't get into table when not analyzed !!!
    
    % loop over paws
    for j = find([kinData.isTrialAnalyzed])
        
        % flip everything relative to first modified paw
        if kinData(j).firstModPaw==referenceModPaw; pawSequence = [1 2 3 4]; else; pawSequence = [4 3 2 1]; end
        flatBin = sessionBins & [flat.trial]==j;
        
        for k = 1:4
            
            % starting paw position
            flat.(['modStepStart_paw' num2str(k)])(flatBin) = kinData(j).modifiedLocations{pawSequence(k)}(1,1,1);
            
            % is in stance
            flat.(['isStance_paw' num2str(k)])(flatBin) = kinData(j).stanceBins(kinData(j).contactInd,pawSequence(k));
            
            % x and z positions
            flat.(['x_paw' num2str(k)])(flatBin) = kinData(j).locations(kinData(j).contactInd, 1, pawSequence(k));
            flat.(['z_paw' num2str(k)])(flatBin) = kinData(j).locations(kinData(j).contactInd, 3, pawSequence(k));
            
            % x and z velocity
            kin = kinData(j).locations(kinData(j).contactInd-velSmps+1:kinData(j).contactInd, :, pawSequence(k));
            flat.(['xVel_paw' num2str(k)])(flatBin) = (kin(end,1)-kin(1,1)) / velSmps / frameRate;
            flat.(['zVel_paw' num2str(k)])(flatBin) = (kin(end,3)-kin(1,3)) / velSmps / frameRate;
        end
    end
end
disp('all done!')






%% PREPARE DATA

% settings
predictors = {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', ...
    'modPawStartingDistance', 'modStepStart', 'isStance', 'x', 'z', 'xVel', 'zVel'};
predictors = {'x_paw2', 'xVel_paw3', 'obsHgt'};
allPaws = false;
outcome = 'isBigStep';


% initializations
cols = flat.Properties.VariableNames;
predictorsTemp = {};
for i = 1:length(predictors)

    % add to list if predictor appears in column as is
    if ismember(predictors{i}, cols)
        predictorsTemp{end+1} = predictors{i};

    % otherwise check if it is a paw specific variable
    elseif ismember([predictors{i} '_paw' num2str(referenceModPaw)], cols)

        if allPaws % add variable for each paw
            for j = 1:4; predictorsTemp{end+1} = [predictors{i} '_paw' num2str(j)]; end        
        else % add only for mod paw
            predictorsTemp{end+1} = [predictors{i} '_paw' num2str(referenceModPaw)];
        end
    end
    
end
predictors = predictorsTemp;



%% set up data
[~, colInds] = ismember(predictors, flat.Properties.VariableNames);
x = table2array(flat(:,colInds));
y = flat.(outcome);

% remove nans
validBins = ~any(isnan([x y]),2);
x = x(validBins,:);
% x = rand(size(x)); % temp
y = y(validBins);

% subsample to even out class sizes
modalClass = mode(y);
inds = find(y==modalClass);
inds = inds(randperm(length(inds), sum(y~=modalClass)));
inds = sort([inds; find(y~=modalClass)]);
x = x(inds,:);
y = y(inds);

xNorm = normalize(x,1); % !!! should not normalize logical vars...

%% MODEL THAT SHIT

% settings
hiddenLayers = [100];



% MULTILAYER PERCEPTRON
net = patternnet(hiddenLayers);
net.divideParam.trainRatio = .7;
net.divideParam.valRatio = .15;
net.divideParam.testRatio = .15;

% train
[net, tr] = train(net, xNorm', y');

% test
outputs = net(xNorm');
errors = gsubtract(y', outputs);
fprintf('\nNEURAL NET train: %.2f, test: %.3f\n', mean(round(outputs(tr.trainInd))==y(tr.trainInd)'), ...
    mean(round(outputs(tr.testInd))==y(tr.testInd)'))
% performance = perform(net, y', outputs)

%% GLM
glm = fitglm(xNorm([tr.trainInd tr.valInd],:), y([tr.trainInd tr.valInd]), ...
    'Distribution', 'binomial', 'VarNames', [predictors outcome]);
fprintf('GLM train: %.2f, test: %.3f\n', mean(round(predict(glm, xNorm(tr.trainInd,:)))==y(tr.trainInd)), ...
    mean(round(predict(glm, xNorm(tr.testInd,:)))==y(tr.testInd)))

%% LASSO
[glmLasso, lassoFit] = lassoglm(xNorm([tr.trainInd tr.valInd],:), y([tr.trainInd tr.valInd]), ...
    'binomial', 'PredictorNames', predictors);
% fprintf('GLM_lasso train: %.2f, test: %.3f\n', mean(round(predict(glmLasso, xNorm(tr.trainInd,:)))==y(tr.trainInd)), ...
%     mean(round(predict(glmLasso, xNorm(tr.testInd,:)))==y(tr.testInd)))

zeroInds = nan(1, length(predictors));
for i = 1:size(glmLasso,1); zeroInds(i) = find(abs(glmLasso(i,:))>0,1,'last'); end
[~, sortInds] = sort(zeroInds, 'descend');
predictorsSorted = predictors(sortInds)

% lassoPlot(glmLasso, lassoFit, 'PlotType', 'Lambda', 'PredictorNames', predictors); legend(predictors)

% !!! check accuracy with different numbers of predictors here...

%% STEPWISE REGRESSION

glmStepwise = stepwiseglm(xNorm([tr.trainInd tr.valInd],:), y([tr.trainInd tr.valInd]), 'constant', ...
    'Upper', 'linear', 'Distribution', 'binomial', 'VarNames', [predictors outcome], 'PEnter', .001);
fprintf('GLM_stepwise train: %.2f, test: %.3f\n', mean(round(predict(glmStepwise, xNorm(tr.trainInd,:)))==y(tr.trainInd)), ...
    mean(round(predict(glmStepwise, xNorm(tr.testInd,:)))==y(tr.testInd)))













