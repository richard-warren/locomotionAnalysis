function [B, t] = makeCosBasis(start, stop, n, varargin)
% create raised cosine basis, (n X t), where each row is a raised cosine
% 'bump', with n centers spaced evenly between start and stop times // the
% time axis is forced to be symmetric around zero (it has odd length and
% the center element is 0)
%
% inspired by: https://github.com/pillowlab/raisedCosineBasis/blob/master/makeRaisedCosBasis.m


% settings
s.matchVal = pi/2;    % (radians) value at which adjacent bases overlap
s.dt = .01;
s.showPlot = false;

% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
centers = linspace(start, stop, n);
dc = centers(2)-centers(1);
f = 2*s.matchVal/dc;  % cosine frequency
tmax = max(abs([start stop])) + pi/f;
t = 0 : s.dt : tmax;
t = [fliplr(-t(2:end)) t];



basisfn = @(t,c) cos(max(min((t-c)*f, pi),-pi)) / 2 + .5;  % time, centers, frequency

B = basisfn(repmat(t,n,1), centers');

if s.showPlot
    figure('color', 'white', 'position', [112.00 967.00 1074.00 173.00]);
    plot(t, B, 'linewidth', 3)
    set(gca, 'box', 'off', 'xlim', [t(1) t(end)])
    xlabel('time (s)')
end