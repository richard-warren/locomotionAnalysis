folder = 'Z:\obstacleData\histology\CER 5\Cer5TiffStack';
addpath('Z:\obstacleData\histology\CER 5\Cer5TiffStack');

%% importing probe traces

tiff_info = imfinfo('MaskProbeTracesLeft.tif');
probeTraceRight = imread('MaskProbeTracesLeft.tif', 1) ; % read in first image

%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('MaskProbeTracesLeft.tif', ii);
    probeTraceRight = cat(3 , probeTraceRight, temp_tiff);
end

%% importing PC layer 

tiff_info = imfinfo('MaskPCLayerLeft.tif');
PCLayerRight = imread('MaskPCLayerLeft.tif', 1) ; % read in first image

%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('MaskPCLayerLeft.tif', ii);
    PCLayerRight = cat(3 , PCLayerRight, temp_tiff);
end

%% importing dentate nucleus

tiff_info = imfinfo('Cer5DentateRight.tif');
dentateRight = imread('Cer5DentateRight.tif', 1) ; % read in first image

%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('Cer5DentateRight.tif', ii);
    dentateRight = cat(3 , dentateRight, temp_tiff);
end

%% test plot in 2d

x = dentateRight(:, :, 13);
test = [];

for i = 1 : size(x, 1)
    for j = 1 : size(x, 2)
        if any(x(i, j))
            test = [test; j, i];
        end
    end
end

figure;
plot(test(:, 1), test(:, 2))
xlim([0, size(x, 2)]);
ylim([0, size(x, 1)]);
set( gca, 'YDir', 'reverse' );

%% downsample the image

% 
% for k = 1:size(dentateRight, 3)
%     for i = 1:size(dentateRight, 2)
%         for j = 1:size(dentateRight, 1)
%             if any(dentateRight(j, i, k))
%                 x = [x; j, i, (k-1)*thickness + 1];
%             end
%         end
%     end
% end
%    

dentate_new = imresize(dentateRight(:, :, 13), 0.3);


for k = 14:28
    temp = imresize(dentateRight(:, :, k), 0.3);
    dentate_new = cat(3 , dentate_new, temp);
end



probe_new = imresize(probeTraceRight(:, :, 13), 0.3);
for k = 14:28
    temp = imresize(probeTraceRight(:, :, k), 0.3);
    probe_new = cat(3, probe_new, temp);
end


PC_new = imresize(PCLayerRight(:, :, 13), 0.3);
for k = 14:28
    temp = imresize(PCLayerRight(:, :, k), 0.3);
    PC_new = cat(3, PC_new, temp);
end


%% test plot in 3d

thickness = 50;
dentateLocation = [];

for k = 1:size(dentate_new, 3)
    for i = 1:size(dentate_new, 2)
        for j = 1:size(dentate_new, 1)
            if any(dentate_new(j, i, k))
                dentateLocation = [dentateLocation; j, i, (k-1)*thickness + 1];
            end
        end
    end
end



thickness = 50;
probeLocation = []

for k = 1:size(probe_new, 3)
    for i = 1:size(probe_new, 2)
        for j = 1:size(probe_new, 1)
            if any(probe_new(j, i, k))
                probeLocation = [probeLocation; i, j, (k-1)*thickness + 1];
            end
        end
    end
end


thickness = 50;
PCLocation = []

for k = 1:size(PC_new, 3)
    for i = 1:size(PC_new, 2)
        for j = 1:size(PC_new, 1)
            if any(PC_new(j, i, k))
                PCLocation = [PCLocation; i, j, (k-1)*thickness + 1];
            end
        end
    end
end


%% plot dentate

s = 20;

scatter3(dentateLocation(:, 2), dentateLocation(:, 1), dentateLocation(:, 3), s, 'filled');
xlim([0, size(dentate_new, 2)]);
ylim([0, size(dentate_new, 1)]);
set( gca, 'YDir', 'reverse' );



hold on;

s = 10;
c = [1, 0.43, 0.54];

scatter3(probeLocation(:, 1), probeLocation(:, 2), probeLocation(:, 3), s, c, 'filled');
xlim([0, size(probe_new, 2)]);
ylim([0, size(probe_new, 1)]);
set( gca, 'YDir', 'reverse' );


xlabel('ML axis')
ylabel('DV axis')
zlabel('AP axis')



hold on;

s = 10;
c = [0.815, 0.56, 0.96];

scatter3(PCLocation(:, 1), PCLocation(:, 2), PCLocation(:, 3), s, c, 'filled');
xlim([0, size(probe_new, 2)]);
ylim([0, size(probe_new, 1)]);
set( gca, 'YDir', 'reverse' );




%% 3D Linear Regression
N = size(probeLocation, 1);                  % number of samples
probeAve = mean(probeLocation,1);            % mean; line of best fit will pass through this point  
dX=bsxfun(@minus,probeLocation,probeAve);    % residuals
C=(dX'*dX)/(N-1);                            % variance-covariance matrix of X
[R,D]=svd(C,0); 

% Coefficient of determineation; R^2 = (explained variance)/(total variance)
D=diag(D);
R2=D(1)/sum(D);

% End-points of a best-fit line (segment); used for visualization only 
x=dX*R(:,1);    % project residuals on R(:,1) 
x_min=min(x);
x_max=max(x);
dx=x_max-x_min;
Xa=(x_min-0.05*dx)*R(:,1)' + probeAve;
Xb=(x_max+0.05*dx)*R(:,1)' + probeAve;
X_end=[Xa;Xb];

%% Plot fitted line

hold on
plot3(X_end(:,1),X_end(: ,2),X_end(:,3),'-g','LineWidth',3) % best fit line 
















