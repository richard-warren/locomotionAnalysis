function checkObsLight(session)

% for a given session shows image montages of the obstacle (as seen in the
% wecam view) on trials when the light is on and off // use to diagnose
% whether the light on the osbtacle is working

try
    load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
        'isLightOn', 'obsOnTimes', 'obsLightOnTimes', 'webCamTimeStamps')

    if any(isLightOn)
        vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'webcam.avi'));
        frames = nan(vid.Height, vid.Height, length(obsOnTimes));

        for j = 1:length(obsOnTimes)
            if isLightOn(j)
                ind = find(webCamTimeStamps>(obsLightOnTimes(knnsearch(obsLightOnTimes, obsOnTimes(j)))), 1, 'first');
            else
                ind = find(webCamTimeStamps>(obsOnTimes(j)), 1, 'first');
            end
            frame = read(vid, ind+2); % +2 to ensure the frame is not missed
            frames(:,:,j) = squeeze(frame(:, end-vid.Height+1:end,1));
        end
        figure('Name', session, 'Position', [37 83 1778 873], 'Color', 'white', 'MenuBar', 'none');
        subplot(1,2,1); montage(frames(:,:,~isLightOn), 'DisplayRange', [0 255]);
        subplot(1,2,2); montage(frames(:,:,isLightOn), 'DisplayRange', [0 255]);
    end
catch
    disp([session ': unable to show obstacle light on vs. off trials'])
end






