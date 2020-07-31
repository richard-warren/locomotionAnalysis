function [accuracies, net, glm] = trainDecisionModels(X, y, iterations, isCategorical, opts)
% train binary classifier using both neural net and GLM, repeating
% iterations times with different train / test splits // returns 2 X
% iterations matrix, where each column is accuracy for NN and GLM fora
% particular train/test split // also returns the neural net and glm from
% only the LAST iteration


s.trainValidationTest = [.7 .15 .15]; % must add to 1
s.hiddenUnits = 100;
s.verbose = false;


% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


accuracies = nan(2, iterations);

for i = 1:iterations

    % NEURAL NETWORK
    net = patternnet(s.hiddenUnits);
    net.divideParam.trainRatio = s.trainValidationTest(1);
    net.divideParam.valRatio = s.trainValidationTest(2);
    net.divideParam.testRatio = s.trainValidationTest(3);
    [net, tr] = train(net, X', y');
    outputs = net(X');
    accuracies(1,i) = mean(round(outputs(tr.testInd))==y(tr.testInd)');
    if s.verbose
        fprintf('\n%i) NEURAL NET train: %.2f, test: %.3f\n', i, mean(round(outputs(tr.trainInd))==y(tr.trainInd)'), ...
            mean(round(outputs(tr.testInd))==y(tr.testInd)'))
    end

    
    % GLM
    glm = fitglm(X([tr.trainInd tr.valInd],:), y([tr.trainInd tr.valInd]), ...
        'Distribution', 'binomial', 'CategoricalVars', isCategorical);
    accuracies(2,i) = mean(round(predict(glm, X(tr.testInd,:)))==y(tr.testInd));
    if s.verbose
        fprintf('%i) GLM train: %.2f, test: %.3f\n', mean(round(predict(glm, X(tr.trainInd,:)))==y(tr.trainInd)), ...
            mean(round(predict(glm, X(tr.testInd,:)))==y(tr.testInd)))
    end
end





