function trainNn(class)

load([getenv('OBSDATADIR') 'svm\trainingData\' class '\labeledFeatures.mat'], 'features', 'labels')
x = features;
I = eye(max(labels));
t = I(labels,:)';

trainFcn = 'trainbr';  % Scaled conjugate gradient backpropagation.

% Create a Pattern Recognition Network
hiddenLayerSize = 100;
net = patternnet(hiddenLayerSize, trainFcn);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y)
tind = vec2ind(t);
yind = vec2ind(y);
percentErrors = sum(tind ~= yind)/numel(tind);

% View the Network
% view(net)

save([getenv('OBSDATADIR') 'svm\classifiers\pawBot2Nn'], 'net')