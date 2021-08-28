% scrape daily logs for baseline sessions, collecting velocity and distance
% traveled.

fname = fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx');
opts = detectImportOptions(fname, 'Sheet', 'dailyLog', 'TextType', 'char');  % need to specify TextType as char; otherwise strange results...
dailyLog = readtable(fname, opts, 'Sheet', 'dailyLog', ...
    'ReadRowNames', true, 'ReadVariableNames', true);
mice = dailyLog.Properties.VariableNames;

%% exclude mice
cohortsToExclude = {'vgt', 'vglut', 'pcp', 'cer8'};
mice = mice(~contains(mice, cohortsToExclude));


%% initialize results table
matInit = nan(length(mice), 30);  % [mouse X days]
data = table(matInit, matInit, nan(length(mice), 1), ...
    'RowNames', mice, 'VariableNames', {'distance', 'velocity', 'days_to_train'});

distancePatterns = {'baseline (', ', training'};  % strings surrounding distance travelled
velocityPatterns = {'speed = ', 'speed ', 'vel '};

for m = 1:length(mice)
    mouseLogs = dailyLog{:, mice{m}};

    % find inds for baseline log entries
    rowInds = find(~cellfun(@isempty, mouseLogs) & ...
                   cellfun(@(x) contains(x, distancePatterns{1}), mouseLogs) & ...
                   cellfun(@(x) contains(x, distancePatterns{2}), mouseLogs));
    
   % get distance travelled and vel for each entry
   if ~isempty(rowInds)
       for idx = 1:min(length(rowInds), size(matInit, 2))
            logEntry = mouseLogs{rowInds(idx)};  % log for the day

            % distance travelled
            startIdx = strfind(logEntry, distancePatterns{1}) + ...
                length(distancePatterns{1});
            endIdx = strfind(logEntry, distancePatterns{2}) - 1;
            distanceStr = logEntry(startIdx:endIdx);
            if distanceStr(end) == 'm'
                distance = str2double(distanceStr(1:end-1));
            else
                % convert from wheel rotations to meters
                distance = str2double(distanceStr);
                distance = (5.4 / 9) * distance;
            end
            data{mice{m}, 'distance'}(idx) = distance;
            
            % velocity
            for pat = velocityPatterns
                startIdx = strfind(logEntry, pat{1}) + length(pat{1});
                if ~isempty(startIdx)
                    endIdx = min(startIdx+3, length(logEntry));
                    velStr = logEntry(startIdx:endIdx);
                    vel = str2double(velStr);
%                     if isnan(vel); keyboard; end
                    data{mice{m}, 'velocity'}(idx) = vel;
                end
            end
       end
       
       % days to train
       daysToTrain = find(data{mice{m}, 'distance'} >= 5.4, 1, 'first');
       if ~isempty(daysToTrain)
           data{mice{m}, 'days_to_train'} = daysToTrain;
       end
       
       % fill in final distance to rewards
       dist = data{mice{m}, 'distance'};
       finalBaselineInd = find(~isnan(dist), 1, 'last') - 1;
       data{mice{m}, 'distance'}(finalBaselineInd+1:end) = dist(finalBaselineInd);
       
       % fill in final vels
       vels = data{mice{m}, 'velocity'};
       finalBaselineInd = find(~isnan(vels), 1, 'last') - 1;
       if finalBaselineInd>0
        data{mice{m}, 'velocity'}(finalBaselineInd+1:end) = vels(finalBaselineInd);
       end
   end
end

%%

close all; figure('color', 'white', 'menubar', 'none', ...
    'Position', [649.00 467.00 390.00 717.00]); hold on
c = [0 0.4470 0.7410];
xlims = [1 12];


distMat = data{:, 'distance'};
velMat = data{:, 'velocity'};
x = 1:size(matInit, 2);

% reward distance
subplot(3,1,1); hold on
% plot(x, distMat', 'color', [0 0 0 .1])     % individual mice

mn = nanmean(distMat, 1);
st = nanstd(distMat, 1);
patch([x fliplr(x)], [(-st+mn) fliplr(st+mn)], c, ...
    'FaceAlpha', .1, 'EdgeColor', 'none')  % error bars
plot(x, mn, 'color', c, 'lineWidth', 3)    % mean

set(gca, 'Box', 'off', 'XLim', xlims, 'XTickLabel', [])
% xlabel('training day')
ylabel('distance per reward (m)')

% velocity
subplot(3,1,2); hold on
% plot(x, velMat', 'color', [0 0 0 .1])     % individual mice

mn = nanmean(velMat(:,2:end), 1);
st = nanstd(velMat(:,2:end), 1);
x = x(2:end);
patch([x fliplr(x)], [(-st+mn) fliplr(st+mn)], c, ...
    'FaceAlpha', .1, 'EdgeColor', 'none')  % error bars
plot(x, mn, 'color', c, 'lineWidth', 3)    % mean

set(gca, 'Box', 'off')
set(gca, 'XLim', xlims)
xlabel('training day')
ylabel('velocity (m)')

% days to train
subplot(3,1,3); hold on
histogram(data{:, 'days_to_train'}, 'EdgeColor', 'none', 'FaceColor', c)
xlabel('days to train')
ylabel('count')

saveas(gcf, 'E:\lab_files\wheel_paper\learning', 'svg')













