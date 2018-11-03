clear all
close all
clc

view(3);
hold on
grid on

z1 = 1; r1 = 1.5;
z2 = 2; r2 = 2.5;
z3 = 3; r3 = 2;

%circles
plotCircle3D([-z1,0,0], [1,-1,0], r1, 'r', 2);
plotCircle3D([-z2,0,0], [1,-1,0], r2, [1, 148/255, 0], 2);
plotCircle3D([-z3,0,0], [1,-1,0], r3, [1, 208/255, 0], 2);
plotCircle3D([z1,0,0], [1,-1,0], r1, 'r', 2);
plotCircle3D([z2,0,0], [1,-1,0], r2, [1, 148/255, 0], 2);
plotCircle3D([z3,0,0], [1,-1,0], r3, [1, 208/255, 0], 2);

% axis
x=-5:0.01:5;
y=zeros(length(x),1);
z=zeros(length(x),1);
plot3(x, y, z, 'color', 'k', 'lineWidth', 1.5);
plot3(5, 0, 0, '>', 'color', 'k', 'lineWidth', 1.5, 'MarkerFaceColor', 'k');
text(5.2, 0, 0, "z");
z=-3:0.01:3;
x=zeros(length(z),1);
y=zeros(length(z),1);
plot3(x, y, z, 'color', 'k', 'lineWidth', 1.5);
plot3(0, 0, 3, '^', 'color', 'k', 'lineWidth', 1.5, 'MarkerFaceColor', 'k');
text(0, 0, 3.4, "y");
y=-3:0.01:3;
x=zeros(length(y),1);
z=zeros(length(y),1);
plot3(x, y, z, 'color', 'k', 'lineWidth', 1.5);
plot3(0, -3, 0, '>', 'color', 'k', 'lineWidth', 1.5, 'MarkerFaceColor', 'k');
text(0, -3.2, 0, "x");

% radius
plotRadius(z1, r1);
plotRadius(z2, r2);
plotRadius(z3, r3);
plotRadius(-z1, r1);
plotRadius(-z2, r2);
plotRadius(-z3, r3);

% origin
text(-0.4, -0.4 , 0, "0", 'color', 'b', 'FontSize', 15);
% zlabels
text(z1, 0, -0.15, "z1", 'color', 'k', 'FontSize', 10);
text(z2, 0, -0.15, "z2", 'color', 'k', 'FontSize', 10);
text(z3, 0, -0.15, "z3", 'color', 'k', 'FontSize', 10);
text(-z1, 0, -0.15, "-z1", 'color', 'k', 'FontSize', 10);
text(-z2, 0, -0.15, "-z2", 'color', 'k', 'FontSize', 10);
text(-z3, 0, -0.15, "-z3", 'color', 'k', 'FontSize', 10);

% currents
legend("I1", "I2", "I3");

set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
set(gca, 'zticklabel', []);
hold off

function plotRadius(xpos, r)
    z=0:0.1:r;
    x=xpos.*ones(length(z),1);
    y=zeros(length(z),1);
    plot3(x, y, z, 'color', 'k', 'lineStyle', '--', 'lineWidth', 1);
end