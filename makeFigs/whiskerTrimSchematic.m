

% settings
faceWid = 1;
faceHgt = 2;
padWid = .4;
padLoc = [.4 1.25];
noseWid = .15;
noseHgt = .1;
wiskAngles = [.2 1] * (pi/2);
wiskLengths = [1 1.5];
eyePos = [.6 .25];
eyeRad = .1;


% pad layout
pad = {'\alpha', '\beta', '\gamma', '\delta', ...
    'A1', 'B1', 'C1', 'D1', 'E1', ...
    'A2', 'B2', 'C2', 'D2', 'E2', ...
    'A3', 'B3', 'C3', 'D3', 'E3', ...
    'A4', 'B4', 'C4', 'D4', 'E4'};
padX = [2 3 4 5, 1 2 3 4 5, 1 2 3 4 5, 1 2 3 4 5, 1 2 3 4 5];
padY = [1 1 1 1, 2 2 2 2 2, 3 3 3 3 3, 4 4 4 4 4, 5 5 5 5 5];
whiskers = pad;
% whiskers = {'\alpha', '\beta', '\gamma'};
bilateral = true;


% initializations
padX = padX - mean(padX);
padY = padY - mean(padY);
bins = ismember(pad, whiskers);
close all;
figure('menubar', 'none', 'color', 'white', 'position', [584 514 584 534]); hold on


% face outline
t = 0:0.01:pi;
x = faceWid*cos(t);
y = faceHgt*sin(t);
patch(x, y, [102, 70, 34]/255, 'EdgeColor', 'black', 'LineWidth', 3, 'facealpha', .75)

% nose, eyes
t = 0:0.01:2*pi;
x = noseWid*cos(t);
y = noseHgt*sin(t) + faceHgt-.5*noseHgt;
patch(x, y, 'black')
rectangle('position', [eyePos(1)-eyeRad, eyePos(2)-eyeRad, 2*eyeRad, 2*eyeRad], 'curvature', [1 1], 'facecolor', 'black', 'edgecolor', 'none');
rectangle('position', [-eyePos(1)-eyeRad, eyePos(2)-eyeRad, 2*eyeRad, 2*eyeRad], 'curvature', [1 1], 'facecolor', 'black', 'edgecolor', 'none');

% whisker pad
scaling = (padWid/range(padX));
x = padX*scaling + padLoc(1);
y = padY*scaling + padLoc(2);
scatter(x, y, 40, [.15 .15 .15], 'filled')
if bilateral; scatter(-x, y, 40, [.15 .15 .15], 'filled'); end

% whiskers
xFrac = (padX-min(padX)) / range(padX);
yFrac = (padY-min(padY)) / range(padY);
for i = find(bins)
    a = interp1([-1 1], wiskAngles, -xFrac(i) + .1*yFrac(i));
    l = interp1([-1 1], wiskLengths, xFrac(i)-yFrac(i));
    plot([x(i) x(i)+cos(a)*l], [y(i) y(i)+sin(a)*l], ...
        'Color', [.15 .15 .15], 'linewidth', 1.5)
    if bilateral
        plot(-[x(i) x(i)+cos(a)*l], [y(i) y(i)+sin(a)*l], ...
            'Color', [.15 .15 .15], 'linewidth', 1.5)
    end
end

% set axis appearance
daspect([1 1 1])