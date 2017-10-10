function [timesInterp, dataInterp] = interpData(times, data, targetFs)

    % interpolate data sampled at unevenly spaced times with targetFs resolution
    % critically, this function ensures the time grid is spaced s.t. it would intercept at 0 if it were extended back in time
    % this is to ensure that different data streams, each with different times of the first measured value, are all sampled at the same grid points
    %
    % input    times: times of data samples
    %          data: sampled data
    %          targetFs: frequency at which to interpolate data
    %
    % output   timesInterp: times of interpolated data points // grid is spaced s.t it would intercept 0 if extended back in time // lowest and highest values are times immediately to the left and right of lowest and highest 'times'
    %          dataInterp:  linearly interpolated data
    

    % interpolate data
    dt = (1/targetFs);
    timesInterp = 0 : dt : (max(times) + dt); % create evenly spaced time vector starting at 0 // note: important to start the grid at 0 s.t. interpolating different datasets at same fs will result in timeStamps at the same points (ie not offset based on time of first event)
    timesInterp = timesInterp(find(timesInterp>min(times))-1:end);

    dataInterp = interp1(times, data, timesInterp, 'linear');
    
end