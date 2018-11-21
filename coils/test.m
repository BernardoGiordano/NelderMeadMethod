clear all
close all
clc

global mu0
mu0 = 4*pi*1e-7; % vacuum permeability [H/m]

% hidden parameters
hI = 5;
hR = 0.8;
hZ = 0.7;
X = [hI hR hZ];

% objective function
z = -6:0.25:6;
fobj = @(x)((1/mean(Bz(X, z)))*(norm(Bz(x, z) - Bz(X, z))*sqrt(inv(length(z)))));
fobj2d = @(x)((1/mean(Bz2d(X, z)))*(norm(Bz2d(x, z) - Bz2d(X, z))*sqrt(inv(length(z)))));
bound1 = @(x)((x(1)+1)^2+(x(2)-1)^2+(x(3)-1)^2 -1.35^2);
bound2 = @(x)(2*x(1)+x(2)-x(3));

% simplex algorythm
bounds = {};
settings = struct('range', 6, 'step', 0.1, 'slices', length(z), 'dimension', 2);
range = struct('Xmin', 0, 'Xmax', 6, 'Ymin', 0, 'Ymax', 1, 'Zmin', 0, 'Zmax', 1);
stop_conditions = struct('maxFlips', 1000, 'minHalving', 1e-5, 'minMargin', 1e-5);
start_conditions = struct('start', [1 0.1], 'length', 1);
obj = NelderMeadMethod(fobj2d, bounds, stop_conditions, start_conditions, settings, range);

% plot ideal minimum
% plot(X(1), X(2), '.', 'color', 'r', 'lineWidth', 4);
plot3(X(1), X(2), X(3), 'x', 'color', 'g', 'lineWidth', 1.5);

function res = Bz(X, z)
    I = X(1);
    R = X(2);
    Z = X(3);
    res = (I*R.^2)./(((Z-z).^2+R.^2).^(3/2)) + (I*R.^2)./(((Z+z).^2+R.^2).^(3/2));
end

function res = Bz2d(X, z)
    res = Bz([X(1) X(2) 0.8], z);
end