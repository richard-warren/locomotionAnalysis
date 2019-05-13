function [sse, yPredicted] = sseval(params, x, y)

% x = 0:(2*pi/100):2*pi*5;
% f = 1;
% p = pi;
% a = 1;
% noise = rand(size(x)) * .5;

% y1 = a*(sin(f*x) + f*x) + noise;
% y2 = a*(sin(f*x+p) + f*x) + noise;
% y = max(y1,y2);


f = params(1);
p = params(2);
a = params(3);

y1 = a*(sin(f*x) + f*x);
y2 = a*(sin(f*x+p) + f*x);
yPredicted = max(y1,y2);

sse = sum((yPredicted - y).^2);