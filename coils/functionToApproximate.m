close all
clear all
clc

% coil radius [m]
R = [0.7 0.8 0.6];
% coil current [A]
I = [3 5 2];
% coil position on the z axis [m]
Z = [0.4 0.7 0.9];

% the actual script
maxZ = max(Z)*2; % evaluate until 2 times the position of the last coil
step = 0.01; % arbitrary
z = -4:0.01:4; % independant variable
bz = symmetricCoilMagneticField(R, I, Z, z);

% plot results
hold on
grid on
% y axis
plot(z, bz, 'lineWidth', 1);
y = [0 max(bz) * 1.25];
x = [0 0];
plot(x, y, 'color', 'k', 'lineWidth', 1);
plot(0, y(2), '^', 'color', 'k', 'lineWidth', 1, 'MarkerFaceColor', 'k')

% title('Magnetic field for N symmetric coils given Z, R and I')
ylabel('$\tilde{B_{z}}$','Interpreter','latex')
xlabel('z')
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);

% function to evaluate magnetic field on the Z axis
function bz = symmetricCoilMagneticField(R, I, Z, z)
% SYMMETRICCOILMAGNETICFIELD Evaluates the magnetic field on the Z axis
% given coil radius, position and current
    % basic value check
    n = length(R);
    if n ~= length(I) || n ~= length(Z)
        error("Arrays must have the same size")
    end

    mu0 = 4*pi*1e-7; % vacuum permeability [H/m]
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