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

% simplex algorythm
obj = Simplex(fobj, [], [], []);
obj.draw();

% plot ideal minimum
plot3(X(1), X(2), X(3), '.', 'color', 'w', 'lineWidth', 4);

function res = Bz(X, z)
    global mu0
    I = X(1);
    R = X(2);
    Z = X(3);
    res = (mu0/2)*(((I*R.^2)./(((Z-z).^2+R.^2).^(3/2))) + ((I*R.^2)./(((Z+z).^2+R.^2).^(3/2))));
end