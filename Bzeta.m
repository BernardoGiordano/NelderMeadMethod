close all 
clear all
clc

mu0 = 4*pi*1e-7; % vacuum permeability [H/m]
% coil radius [m]
R = [0.7 0.8 0.6];
% coil current [A]
I = [0.3 0.5 0.2];
% coil position on the z axis [m]
Z = [0.4 0.9 1.4];

I2 = I(2);
R2 = R(2);
Z2 = Z(2);
 % evaluate until 2 times the position of the last coil
step = 0.01; % arbitrary
z = -4:step:4; % independant variable
% F = (mu0/2)*(((I2*(R2^2))/((((Z2-z)^2)+(R2^2))^(3/2))) + ((I2*(R2^2))/((((Z2+z)^2)+(R2^2))^(3/2))));
Bz = funzione(R2,I2,Z2,z);


function Bz=funzione(R2, I2, Z2, z)
mu0 = 4*pi*1e-7; % vacuum permeability [H/m]
 Bz = (mu0/2)*(((I2*(R2.^2))./((((Z2-z).^2)+(R2.^2)).^(3/2))) + ((I2*(R2.^2))./((((Z2+z).^2)+(R2.^2)).^(3/2))));
end

    
