close all
clear all
clc

testI = [4.953 4.981 4.954 4.968 4.093 3.406];
testR = [0.791 0.796 0.789 0.794 0.991 0.885];
testZ = [0.699 0.699 0.698 0.699 0.495 0.442];

% coil radius [m]
R = [0.7 0.8 0.6];
% coil current [A]
I = [3 5 2];
% coil position on the z axis [m]
Z = [0.4 0.7 0.9];

% the actual script
maxZ = 3; % evaluate until 2 times the position of the last coil
step = 0.01; % arbitrary
z = -maxZ:0.01:maxZ; % independant variable
bz = symmetricCoilMagneticField(R, I, Z, z);
bz1 = symmetricCoilMagneticField([R(1) testR(1) R(3)], [I(1) testI(1) I(3)], [Z(1) testZ(1) Z(3)], z);
bz2 = symmetricCoilMagneticField([R(1) testR(2) R(3)], [I(1) testI(2) I(3)], [Z(1) testZ(2) Z(3)], z);
bz3 = symmetricCoilMagneticField([R(1) testR(3) R(3)], [I(1) testI(3) I(3)], [Z(1) testZ(3) Z(3)], z);
bz4 = symmetricCoilMagneticField([R(1) testR(4) R(3)], [I(1) testI(4) I(3)], [Z(1) testZ(4) Z(3)], z);
bz5 = symmetricCoilMagneticField([R(1) testR(5) R(3)], [I(1) testI(5) I(3)], [Z(1) testZ(5) Z(3)], z);
bz6 = symmetricCoilMagneticField([R(1) testR(6) R(3)], [I(1) testI(6) I(3)], [Z(1) testZ(6) Z(3)], z);

% plot results
hold on
grid on
% y axis
plot(z, bz, '--', 'lineWidth', 0.5);
plot(z, bz1, 'lineWidth', 1);
plot(z, bz2, 'lineWidth', 1);
plot(z, bz3, 'lineWidth', 1);
plot(z, bz4, 'lineWidth', 1);
plot(z, bz5, 'lineWidth', 1);
plot(z, bz6, 'lineWidth', 1);
legend('Ideale', 'Test 1', 'Test 2', 'Test 3', 'Test 4', 'Test 5', 'Test 6')
%y = [0 max(bz) * 1.25];
%x = [0 0];
%plot(x, y, 'color', 'k', 'lineWidth', 1);
%plot(0, y(2), '^', 'color', 'k', 'lineWidth', 1, 'MarkerFaceColor', 'k')

% title('Magnetic field for N symmetric coils given Z, R and I')
ylabel('$\tilde{B_{z}}$ [T]', 'Interpreter', 'latex', 'fontweight', 'bold', 'fontsize', 13)
xlabel('z [m]')
% set(gca, 'xticklabel', []);
% set(gca, 'yticklabel', []);

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