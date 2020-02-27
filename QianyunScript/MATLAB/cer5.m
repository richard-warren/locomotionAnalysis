%% Cer5 reconstructed data

entryPoint_left1L = [411.6875, 145.5379, 160.305]; entryPoint_left1L = distanceConverter(entryPoint_left1L);
PCPoint_left1L = [385.8273, 115.9455, 428.4403]; PCPoint_left1L = distanceConverter(PCPoint_left1L);

GCChannel_left1L = [378.4557, 107.5101, 504.8737; 377.74, 106.69, 512.29; 374.5194, 103, 545.6877; 373.4459, 101.7772, 556.8187;...
                    372.7302, 100.9583, 564.2395];
GCChannel_left1L = distanceConverter(GCChannel_left1L);

hold on
plot3(entryPoint_left1L(:, 1), entryPoint_left1L(:, 2), entryPoint_left1L(:, 3), '.k', 'MarkerSize', 30);
plot3(PCPoint_left1L(:, 1), PCPoint_left1L(:, 2), PCPoint_left1L(:, 3), '.r', 'MarkerSize', 30);
plot3(GCChannel_left1L(:, 1), GCChannel_left1L(:, 2), GCChannel_left1L(:, 3), '.c', 'MarkerSize', 30);
points = [GCChannel_left1L(5, :); entryPoint_left1L];
plot3(points(:,1),points(:,2),points(:,3),'-g','LineWidth',3)



entryPoint_left2L = [528.3, 167.9345, 139.21]; entryPoint_left2L = distanceConverter(entryPoint_left2L);
GCChannel_left2L = [486.4596, 120.0594, 478.1476; 481.46, 119.59, 486.05]; GCChannel_left2L = distanceConverter(GCChannel_left2L);
hold on
plot3(entryPoint_left2L(:, 1), entryPoint_left2L(:, 2), entryPoint_left2L(:, 3), '.k', 'MarkerSize', 30);
plot3(GCChannel_left2L(:, 1), GCChannel_left2L(:, 2), GCChannel_left2L(:, 3), '.c', 'MarkerSize', 30);
points = [GCChannel_left2L(2, :); entryPoint_left2L];
plot3(points(:,1),points(:,2),points(:,3),'-g','LineWidth',3)



entryPoint_left3L = [539.96, 190.3245, 135.3059]; entryPoint_left3L = distanceConverter(entryPoint_left3L);
GCChannel_left3L = [503.2703, 146.1998, 477.0663; 502.913, 145.77, 480.3947]; GCChannel_left3L = distanceConverter(GCChannel_left3L);
hold on
plot3(entryPoint_left3L(:, 1), entryPoint_left3L(:, 2), entryPoint_left3L(:, 3), '.k', 'MarkerSize', 30);
plot3(GCChannel_left3L(:, 1), GCChannel_left3L(:, 2), GCChannel_left3L(:, 3), '.c', 'MarkerSize', 30);
points = [GCChannel_left3L(2, :); entryPoint_left3L];
plot3(points(:,1),points(:,2),points(:,3),'-g','LineWidth',3)

entryPoint_right1L = [983.7, 155.6, 168]; entryPoint_right1L = distanceConverter(entryPoint_right1L);
PCPoint_right1L = [952.2, 123.8, 437.6]; PCPoint_right1L = distanceConverter(PCPoint_right1L);
GCChannel_right1L = [949.9, 121.4, 457.6; 945.2403, 116.6979, 497.5162]; GCChannel_right1L = distanceConverter(GCChannel_right1L);
hold on
plot3(entryPoint_right1L(:, 1), entryPoint_right1L(:, 2), entryPoint_right1L(:, 3), '.k', 'MarkerSize', 30);
plot3(PCPoint_right1L(:, 1), PCPoint_right1L(:, 2), PCPoint_right1L(:, 3), '.r', 'MarkerSize', 30);

plot3(GCChannel_right1L(:, 1), GCChannel_right1L(:, 2), GCChannel_right1L(:, 3), '.g', 'MarkerSize', 30);
points = [GCChannel_right1L(2, :); entryPoint_right1L];
plot3(points(:,1),points(:,2),points(:,3),'-c','LineWidth',3)



entryPoint_right1R = [1019.6973, 156.1673, 183]; entryPoint_right1R = distanceConverter(entryPoint_right1R);
PCPoint_right1R = [998.0118, 134.262, 368.7273; 988.2173, 124.3692, 452.6042]; PCPoint_right1R = distanceConverter(PCPoint_right1R);
GCChannel_right1R = [980.0571, 116.1252, 522.5016; 979.6685, 115.7327, 525.83; 978.8912, 114.9475, 532.4869; ...
                     978.114, 114.1624, 539.1438; 977.7253, 113.7698, 542.4723]; 
GCChannel_right1R = distanceConverter(GCChannel_right1R);

hold on
plot3(entryPoint_right1R(:, 1), entryPoint_right1R(:, 2), entryPoint_right1R(:, 3), '.k', 'MarkerSize', 30);
plot3(PCPoint_right1R(:, 1), PCPoint_right1R(:, 2), PCPoint_right1R(:, 3), '.r', 'MarkerSize', 30);
plot3(GCChannel_right1R(:, 1), GCChannel_right1R(:, 2), GCChannel_right1R(:, 3), '.g', 'MarkerSize', 30);
points = [GCChannel_right1R(2, :); entryPoint_right1R];
plot3(points(:,1),points(:,2),points(:,3),'-c','LineWidth',3)



entryPoint_right2L = [1070.2595, 170.4, 250]; entryPoint_right2L = distanceConverter(entryPoint_right2L);
PCPoint_right2L = [1062.7357, 162.8, 314.4387; 1050.0664, 150, 422.9461]; PCPoint_right2L = distanceConverter(PCPoint_right2L);
GCChannel_right2L = [1040.7237, 140.5624, 502.962]; GCChannel_right2L = distanceConverter(GCChannel_right2L);

hold on
plot3(entryPoint_right2L(:, 1), entryPoint_right2L(:, 2), entryPoint_right2L(:, 3), '.k', 'MarkerSize', 30);
plot3(PCPoint_right2L(:, 1), PCPoint_right2L(:, 2), PCPoint_right2L(:, 3), '.r', 'MarkerSize', 30);
plot3(GCChannel_right2L(:, 1), GCChannel_right2L(:, 2), GCChannel_right2L(:, 3), '.g', 'MarkerSize', 30);
points = [GCChannel_right2L(1, :); entryPoint_right2L];
plot3(points(:,1),points(:,2),points(:,3),'-c','LineWidth',3)



entryPoint_right2R = [1106.2423, 171, 275]; entryPoint_right2R = distanceConverter(entryPoint_right2R);
PCPoint_right2R = [1087.9767, 152.557, 431.437]; PCPoint_right2R = distanceConverter(PCPoint_right2R);
GCChannel_right2R = [1082.5359, 147.061, 478.0353; 1082.1472, 146.6684, 481.3637]; GCChannel_right2R = distanceConverter(GCChannel_right2R);
hold on
plot3(entryPoint_right2R(:, 1), entryPoint_right2R(:, 2), entryPoint_right2R(:, 3), '.k', 'MarkerSize', 30);
plot3(PCPoint_right2R(:, 1), PCPoint_right2R(:, 2), PCPoint_right2R(:, 3), '.r', 'MarkerSize', 30);
plot3(GCChannel_right2R(:, 1), GCChannel_right2R(:, 2), GCChannel_right2R(:, 3), '.g', 'MarkerSize', 30);
points = [GCChannel_right2R(2, :); entryPoint_right2R];
plot3(points(:,1),points(:,2),points(:,3),'-c','LineWidth',3)


entryPoint_right3R = [1078.1506, 180.192, 242.1994]; entryPoint_right3R = distanceConverter(entryPoint_right3R);
hold on
plot3(entryPoint_right3R(:, 1), entryPoint_right3R(:, 2), entryPoint_right3R(:, 3), '.k', 'MarkerSize', 30);


%%


mouseID = 'cer5';


ephysHistoData = cell(8, 6);

ephysHistoData{1, 1} ='191007_003';
ephysHistoData{1, 2} = mouseID;
ephysHistoData{1, 3} = 'right';
ephysHistoData{1, 4} = entryPoint_right1L;
ephysHistoData{1, 5} = PCPoint_right1L;
ephysHistoData{1, 6} = GCChannel_right1L;

ephysHistoData{2, 1} ='191007_003';
ephysHistoData{2, 2} = mouseID;
ephysHistoData{2, 3} = 'right';
ephysHistoData{2, 4} = entryPoint_right1R;
ephysHistoData{2, 5} = PCPoint_right1R;
ephysHistoData{2, 6} = GCChannel_right1R;


ephysHistoData{3, 1} = '191008_003';
ephysHistoData{3, 2} = mouseID;
ephysHistoData{3, 3} = 'right';
ephysHistoData{3, 4} = entryPoint_right2L;
ephysHistoData{3, 5} = PCPoint_right2L;
ephysHistoData{3, 6} = GCChannel_right2L;


ephysHistoData{4, 1} = '191008_003';
ephysHistoData{4, 2} = mouseID;
ephysHistoData{4, 3} = 'right';
ephysHistoData{4, 4} = entryPoint_right2R;
ephysHistoData{4, 5} = PCPoint_right2R;
ephysHistoData{4, 6} = GCChannel_right2R;

ephysHistoData{5, 1} = '191009_003';
ephysHistoData{5, 2} = mouseID;
ephysHistoData{5, 3} = 'left';
ephysHistoData{5, 4} = entryPoint_left1L;
ephysHistoData{5, 5} = PCPoint_left1L;
ephysHistoData{5, 6} = GCChannel_left1L;

ephysHistoData{6, 1} = '191009_003';
ephysHistoData{6, 2} = mouseID;
ephysHistoData{6, 3} = 'left';
ephysHistoData{6, 4} = [];
ephysHistoData{6, 5} = [];
ephysHistoData{6, 6} = [];


ephysHistoData{7, 1} = '191010_003';
ephysHistoData{7, 2} = mouseID;
ephysHistoData{7, 3} = 'left';
ephysHistoData{7, 4} = entryPoint_left2L;
ephysHistoData{7, 5} = [];
ephysHistoData{7, 6} = GCChannel_left2L;

ephysHistoData{8, 1} = '191011_000';
ephysHistoData{8, 2} = mouseID;
ephysHistoData{8, 3} = 'left';
ephysHistoData{8, 4} = entryPoint_left3L;
ephysHistoData{8, 5} = [];
ephysHistoData{8, 6} = GCChannel_left3L;

%%
fileName = fullfile('D:\WorkFromHome_02212020\ephysHistoData', 'ephysHistoData.mat');
if ~exist(fileName)
    save(fileName, 'ephysHistoData');
else
    temp = load(fileName);
    ephysHistoDataOld = temp.ephysHistoData;
    sessionsOld = ephysHistoDataOld(:, 1);
    sessions = ephysHistoData(:, 1);
    inds = zeros(length(sessionsOld), length(sessions));
    for i = 1:length(sessions)
        inds(:, i) = strcmp(sessions{i}, sessionsOld);
    end
    
    if sum(sum(inds)) == 0
        ephysHistoData = [ephysHistoDataOld; ephysHistoData];
        save(fileName, 'ephysHistoData');
    else
        [updateIndsOld, updateIndsNew] = find(inds == 1);
        for i = 1:length(updateIndsOld)
            ephysHistoDataOld(updateIndsOld(i), :) = ephysHistoData(updateIndsNew(i), :);
            ephysHistoData(updateIndsNew(i), :) = [];
        end
        ephysHistoData = [ephysHistoDataOld; ephysHistoData];
        save(fileName, 'ephysHistoData');
    end    
end






















