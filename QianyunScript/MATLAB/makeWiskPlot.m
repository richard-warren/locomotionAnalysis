
%% create whisker pattern (with correct color pattern)

x = 1:7;
y = 1:3;
[x_value, y_value] = meshgrid(x, y);
sz = 200;

c = jet(9);

figure('Color', 'white');


for i = 1:3
    y_value = repmat(i, 1, length(x_value));
    scatter(x_value(i, :), y_value, sz, c((2*i-1), :), 'filled');
    hold on;
end

hold on;
scatter(8, 1.5, sz, c(2,:), 'filled');
hold on;
scatter(8, 2.5, sz, c(4,:), 'filled');
hold on;
scatter(8, 3.5, sz, c(6,:), 'filled');
hold on;
scatter(8, 4.5, sz, c(8,:), 'filled');
set(gca, 'visible', 'off');


x = 4:7;
y = 4:5;
[x_value, y_value] = meshgrid(x, y);
hold on

for i = 1:2
    y_value = repmat(i+3, 1, length(x_value));
    scatter(x_value(i, :), y_value, sz, c((2*i+5), :), 'filled');
    hold on;
end

set(gca, 'visible', 'off');


%% plot the whiskers (straight line version)

x_value = [1:8; 1:8; 1:8];

y_value = [1 1 1 1 1 1 1 1.5; 2 2 2 2 2 2 2 2.5; 3 3 3 3 3 3 3 3.5];

% wiskend_x = x_value - 3;
% wiskend_y = y_value - 3;

% hold on
% plot([x_value(:)'; wiskend_x(:)'], [y_value(:)'; wiskend_y(:)'], '-', 'Color', [0.5 0.5 0.5 0.5])

for i = 1:3
    
    wiskend_x = x_value - (3.5-i);
    wiskend_y = y_value - (3-0.3*i);
    
    x = [x_value(i, 1:7); wiskend_x(i, 1:7)];
    y = [y_value(i, 1:7); wiskend_y(i, 1:7)];
    hold on
    plot(x, y, '-', 'Color',  c((2*i-1), :))
end



for i = 1:3
    
    wiskend_x = x_value - (4-i);
    wiskend_y = y_value - (3+0.3*i);
    
    x = [x_value(i, 8); wiskend_x(i, 8)];
    y = [y_value(i, 8); wiskend_y(i, 8)];
    hold on
    plot(x, y, '-', 'Color',  c((2*i), :))
end



%% %% plot the whiskers (curved line version)



x_value = [1:8; 1:8; 1:8];

y_value = [1 1 1 1 1 1 1 1.5; 2 2 2 2 2 2 2 2.5; 3 3 3 3 3 3 3 3.5];

r = 7;
n = 25

circr = @(radius, rad_ang, p0)  [radius*cos(rad_ang)+p0(1);  radius*sin(rad_ang)+p0(2)];         % Circle Function For Angles In Radians


for i = 1:3
    angle = linspace((9+2*i)*pi/(5+i) , 2*pi, n);
    
    for j = 1:7
        p0 = [x_value(i, j)-r, y_value(i, j)];
        wiskend_xy = circr(r, angle, p0);
        
        hold on
        plot(wiskend_xy(1, :), wiskend_xy(2,:), '-', 'Color',  c((2*i-1), :))
    end
end


for i = 1:3
    angle = linspace((7+2*i)*pi/(4+i) , 2*pi, n);
    
    
    p0 = [x_value(i, 8)-r, y_value(i, 8)];
    wiskend_xy = circr(r, angle, p0);
    
    hold on
    plot(wiskend_xy(1, :), wiskend_xy(2,:), '-', 'Color',  c((2*i), :))
    
end



%% generate the other whisker pad

x = 1:7;
y = 1:3;
[x_value, y_value] = meshgrid(x, y);

c = parula(8);

x_value = -x_value;
y_value = -y_value;

%hold on;
figure('Color', 'white');

for i = 1:3
    y_value = repmat(i, 1, length(x_value));
    scatter(x_value(i, :), y_value, [], c((2*i-1), :), 'filled');
    hold on;
end

hold on;
scatter(-8, 1.5, [], c(2,:), 'filled');
hold on;
scatter(-8, 2.5, [], c(4,:), 'filled');
hold on;
scatter(-8, 3.5, [], c(6,:), 'filled');
set(gca, 'visible', 'off');


x_value = [1:8; 1:8; 1:8];

y_value = [1 1 1 1 1 1 1 1.5; 2 2 2 2 2 2 2 2.5; 3 3 3 3 3 3 3 3.5];

r = 7;
n = 25

circr = @(radius, rad_ang, p0)  [radius*cos(rad_ang)+p0(1);  radius*sin(rad_ang)+p0(2)];         % Circle Function For Angles In Radians


for i = 1:3
    angle = linspace((9+2*i)*pi/(5+i) , 2*pi, n);
    
    for j = 1:7
        p0 = [x_value(i, j)-r, y_value(i, j)];
        wiskend_xy = circr(r, angle, p0);
        
        hold on
        plot(-wiskend_xy(1, :), wiskend_xy(2,:), '-', 'Color',  c((2*i-1), :))
    end
end


for i = 1:3
    angle = linspace((7+2*i)*pi/(4+i) , 2*pi, n);
    
    
    p0 = [x_value(i, 8)-r, y_value(i, 8)];
    wiskend_xy = circr(r, angle, p0);
    
    hold on
    plot(-wiskend_xy(1, :), wiskend_xy(2,:), '-', 'Color',  c((2*i), :))
    
end
























