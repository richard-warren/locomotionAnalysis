% Addpath of the often use folders
setenv('GITDIR', 'D:\DCN_project\Github');

addpath(genpath('D:\DCN_project\Github'));
addpath('Z:\obstacleData\ephys');
addpath('Z:\obstacleData\spreadSheets');

warning('off','all');
rmpath(genpath('D:\DCN_project\Github\LocomotionAnalysis_QZ'));
warning('on','all');


display('finish loading startup.m');