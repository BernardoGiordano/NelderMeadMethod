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


z = -2:0.01:2;
fobj = @(x)(norm(Bz(x, z) - Bz(X, z), 2));
obj = Simplex(fobj);
plot3(0.5, 0.8, 0.9, '^', 'color', 'w');
obj.draw();

function res = Bz(V, z)
    global mu0
    I = V(1);
    R = V(2);
    Z = V(3);
    res = (mu0/2)*(((I*R.^2)./(((Z-z).^2+R.^2).^(3/2))) + ((I*R.^2)./(((Z+z).^2+R.^2).^(3/2))));
end