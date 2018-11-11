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
z = -4:0.01:4;
fobj = @(x)(norm(Bz(x, z) - Bz(X, z), 2));
bound1 = @(x)((x(1)+1)^2+(x(2)-1)^2+(x(3)-1)^2 -1.35^2);
bound2 = @(x)(2*x(1)+x(2)-x(3));

% simplex algorythm
bounds = {bound1, bound2};
stop_conditions = struct('steps', 1000, 'minArea', 1e-5, 'minHalving', 1e-5, 'minMargin', 1e-5);
start_conditions = struct('start', [0.0, 0.0, 0.0], 'area', 0.5);
obj = NelderMeadMethod(fobj, bounds, stop_conditions, start_conditions);

% plot ideal minimum
plot3(X(1), X(2), X(3), '.', 'color', 'w', 'lineWidth', 4);
plot3(0, 0, 0, '.', 'color', 'w', 'lineWidth', 8);

function res = Bz(X, z)
    global mu0
    I = X(1);
    R = X(2);
    Z = X(3);
    res = (mu0/2)*(((I*R.^2)./(((Z-z).^2+R.^2).^(3/2))) + ((I*R.^2)./(((Z+z).^2+R.^2).^(3/2))));
end