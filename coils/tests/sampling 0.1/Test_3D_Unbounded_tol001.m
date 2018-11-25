addpath(genpath('../'));

clear all
close all
clc

% hidden parameters
hI = 5;
hR = 0.8;
hZ = 0.7;
X = [hI hR hZ];

%%%%%%%%%%%%%% EDIT ME %%%%%%%%%%%%%%
test_params = struct( ...
    'dimension', 3, ...
    'plot', false, ...
    'minimum', X, ...
    'sampling_step', 0.1, ...
    'sampling_range', 3, ...
    'sample_amount', -1, ...
    'max_flips', 10000, ...
    'tolerance', 0.01, ...
    'minLength', 1e-3, ...
    'start_point', [3 0.3 0.5], ...
    'length', 0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = -test_params.sampling_range:test_params.sampling_step:test_params.sampling_range;
test_params.sample_amount = length(z);
disp("Test parameters")
disp(test_params)

%%%%%%%%%%%%%% EDIT ME %%%%%%%%%%%%%%
bounds = {};
fobj = @(x)((1/mean(Bz(X, z)))*(norm(Bz(x, z) - Bz(X, z))*sqrt(inv(length(z)))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simplex algorythm
settings = struct('step', test_params.sampling_step, 'slices', floor(length(z)/2), 'plot', test_params.plot, 'dimension', test_params.dimension);
range = struct('Xmin', 0, 'Xmax', 6, 'Ymin', 0, 'Ymax', 1, 'Zmin', 0, 'Zmax', 1);
stop_conditions = struct('maxFlips', test_params.max_flips, 'tolerance', test_params.tolerance, 'minLength', test_params.minLength);
start_conditions = struct('start', test_params.start_point, 'length', test_params.length);
obj = NelderMeadMethod(fobj, bounds, stop_conditions, start_conditions, settings, range);
% label
xlabel('I [A]');
ylabel('R [m]');
zlabel('Z [m]');
% display results
disp("Results")
disp(obj.getResults());
if test_params.plot
    % plot ideal minimum
    plot3(X(1), X(2), X(3), 'x', 'color', 'y', 'lineWidth', 1.5);
end