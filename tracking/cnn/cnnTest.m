


% load alexNet
net = alexnet;

%% rick mouse
vid = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\180109_002\runTop.mp4');
img = read(vid, 10000);

%% mouse
img = imread('C:\Users\rick\Desktop\mouse.jpg');

%% cat
img = imread('C:\Users\rick\Desktop\cat.jpg');

%% classify
imgResized = imresize(img, [227 227]);
label = classify(net, imgResized)