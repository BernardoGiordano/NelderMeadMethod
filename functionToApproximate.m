close all
clear all
clc

% coil radius [m]
R = [1 1 1];
% coil current [A]
I = [1 1 1];
% coil position on the z axis [m]
Z = [1 1 1];

% the actual script
maxZ = max(Z)*4; % evaluate until 4 times the position of the last coil
step = 0.01; % arbitrary
z = -maxZ:step:maxZ; % independant variable
bz = symmetricCoilMagneticField(R, I, Z, z);
% plot results
hold on
grid on
plot(z, bz);
title('Magnetic field for N symmetric coils given Z, R and I')
ylabel('Bz')
xlabel('z')


% function to evaluate magnetic field on the Z axis
function bz = symmetricCoilMagneticField(R, I, Z, z)
% SYMMETRICCOILMAGNETICFIELD Evaluates the magnetic field on the Z axis
% given coil radius, position and current
    % basic value check
    n = length(R);
    if n ~= length(I) || n ~= length(Z)
        error("Arrays must have the same size")
    end

    mu0 = 4*pi*10e-7; % vacuum permeability [H/m]
    iterations = length(z); % how many iteration we need
    bz = zeros(iterations, 1); % output array, preallocating memory for efficiency
    
    for k = 1:iterations
        inc = 0;
        dz = z(k);
        for i = 1:n
            inc = inc + (R(i).^2 * I(i).^2)./((Z(i) - dz).^2 + R(i).^2).^(3/2);
            % also add symmetric portion in the same cycle
            inc = inc + (R(i).^2 * I(i).^2)./((Z(i) + dz).^2 + R(i).^2).^(3/2);
        end
        bz(k) = inc;
    end
    bz = mu0 / 2 * bz;
end