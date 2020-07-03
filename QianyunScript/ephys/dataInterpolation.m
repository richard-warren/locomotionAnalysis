function [data1Interp, time1Interp, data2, time2] = dataInterpolation(data1, time1, data2, time2, opts)


% settings
s.data1Freq = 250; % Hz
s.data2Freq = 1000; % Hz
s.interpMethod = 'spline'; % interpolation method

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end



if time1(1) <= time2(1); startTime = time2(1); else; startTime = time1(1); end;
if time1(end) >= time2(end); endTime = time2(end); else; endTime = time1(end); end;

if size(time1, 2) ~= 1; time1 = time1'; end
if size(time2, 2) ~= 1; time2 = time2'; end


ind(1) = knnsearch(time2, startTime);
ind(2) = knnsearch(time2, endTime);
time2 = time2(ind(1):ind(2));
data2 = data2(ind(1):ind(2));

clear ind
ind(1) = knnsearch(time1,startTime);
ind(2) = knnsearch(time1, endTime);
time1 = time1(ind(1):ind(2));
data1 = data1(ind(1):ind(2));

time1Interp = linspace(time1(1), time1(end), length(data2));
time1Interp = time1Interp';
indsInterp = linspace(ind(1), ind(2), length(time1Interp));

if islogical(data1)
    data1 = double(data1);
end

if size(data1, 2) == 1; data1 = data1'; end
data1Interp = interp1(1:length(data1), data1, indsInterp, s.interpMethod);


if length(data1Interp) ~= length(data2)
    warning('Data length DOES NOT match!!!! -> length(data1Interp) ~= length(data2)')
end  

if mean(abs(time1Interp - time2)) > 0.001
    warning('Time match FAILED somehow...-> mean(abs(time1Interp - time2)) > 0.0001');
end


pause(.001)
end
