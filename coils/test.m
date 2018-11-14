clear all
close all
clc

global mu0
mu0 = 4*pi*1e-7; % vacuum permeability [H/m]

% hidden parameters
hI = 0.5;
hR = 0.8;
hZ = 0.9;
X = [hI hR hZ];

% objective function
z = -1:0.1:1;
Bzh = mean(Bz(X, z));
fobj = @(x)((1/Bzh)*(norm(Bz(x, z) - Bz(X, z)).*sqrt(inv(length(z)))));
fobj2d = @(x)(norm(Bz2d(x, z) - Bz2d(X, z)));
bound1 = @(x)((x(1)+1)^2+(x(2)-1)^2+(x(3)-1)^2 -1.35^2);
bound2 = @(x)(2*x(1)+x(2)-x(3));

% simplex algorythm
bounds = {};
settings = struct('range', 1, 'step',0.01, 'slices', 20, 'func_counter', 0, 'dimension', 3);
stop_conditions = struct('maxFlips', 1000, 'minHalving', 1e-5, 'minMargin', 1e-5);
start_conditions = struct('start', [0.4 0.1 0.5], 'length', 0.5);
obj = NelderMeadMethod(fobj, bounds, stop_conditions, start_conditions, settings);

% plot ideal minimum
% plot(X(1), X(2), '.', 'color', 'r', 'lineWidth', 4);
% plot3(X(1), X(2), X(3), '.', 'color', 'r', 'lineWidth', 4);

function res = Bz(X, z)
    global mu0
    I = X(1);
    R = X(2);
    Z = X(3);
    res = (mu0/2)*(((I*R.^2)./(((Z-z).^2+R.^2).^(3/2))) + ((I*R.^2)./(((Z+z).^2+R.^2).^(3/2))));
end

function res = Bz2d(X, z)
    global mu0
    I = X(1);
    R = X(2);
    Z = 0.8;
    res = (mu0/2)*(((I*R.^2)./(((Z-z).^2+R.^2).^(3/2))) + ((I*R.^2)./(((Z+z).^2+R.^2).^(3/2))));
end