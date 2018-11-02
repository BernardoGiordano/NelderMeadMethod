clear all
close all
clc

view(3);

hold on
grid on
xlabel("z");
ylabel("x");
zlabel("y");

%circles
plotCircle3D([-1,0,0], [1,-1,0], 1.5, 'r', 2);
plotCircle3D([-2,0,0], [1,-1,0], 2.5, [1, 148/255, 0], 2);
plotCircle3D([-3,0,0], [1,-1,0], 2, [1, 208/255, 0], 2);
plotCircle3D([1,0,0], [1,-1,0], 1.5, 'r', 2);
plotCircle3D([2,0,0], [1,-1,0], 2.5, [1, 148/255, 0], 2);
plotCircle3D([3,0,0], [1,-1,0], 2, [1, 208/255, 0], 2);

% axis
x=-5:0.01:5;
y=zeros(length(x),1);
z=zeros(length(x),1);
plot3(x, y, z, 'color', 'k', 'lineWidth', 1.5);
plot3(5, 0, 0, '>', 'color', 'k', 'lineWidth', 1.5, 'MarkerFaceColor', 'k')
z=-3:0.01:3;
x=zeros(length(z),1);
y=zeros(length(z),1);
plot3(x, y, z, 'color', 'k', 'lineWidth', 1.5);
plot3(0, 0, 3, '^', 'color', 'k', 'lineWidth', 1.5, 'MarkerFaceColor', 'k')
y=-3:0.01:3;
x=zeros(length(y),1);
z=zeros(length(y),1);
plot3(x, y, z, 'color', 'k', 'lineWidth', 1.5);
plot3(0, -3, 0, '>', 'color', 'k', 'lineWidth', 1.5, 'MarkerFaceColor', 'k')


% radius
plotRadius(1, 1.5);
plotRadius(2, 2.5);
plotRadius(3, 2);
plotRadius(-1, 1.5);
plotRadius(-2, 2.5);
plotRadius(-3, 2);

text(-0.4, -0.4 , 0, "0", 'color', 'b', 'FontSize', 15);

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