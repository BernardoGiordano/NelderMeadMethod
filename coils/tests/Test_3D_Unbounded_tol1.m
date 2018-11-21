clear all
close all
clc

% hidden parameters
hI = 5;
hR = 0.8;
hZ = 0.7;
X = [hI hR hZ];

disp("Test parameters")
%%%%%%%%%%%%%% EDIT ME %%%%%%%%%%%%%%
test_params = struct( ...
    'dimension', 3, ...
    'minimum', X, ...
    'sampling_step', 0.1, ...
    'sampling_range', 2, ...
    'sample_amount', -1, ...
    'max_flips', 1000, ...
    'tolerance', 1, ...
    'start_point', [3 0.3 0.5], ...
    'length', 0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = -test_params.sampling_range:test_params.sampling_step:test_params.sampling_range;
test_params.sample_amount = length(z);
disp(test_params)

%%%%%%%%%%%%%% EDIT ME %%%%%%%%%%%%%%
bounds = {};
fobj = @(x)((1/mean(Bz(X, z)))*(norm(Bz(x, z) - Bz(X, z))*sqrt(inv(length(z)))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simplex algorythm
settings = struct('step', test_params.sampling_step, 'slices', length(z), 'dimension', test_params.dimension);
range = struct('Xmin', 0, 'Xmax', 6, 'Ymin', 0, 'Ymax', 1, 'Zmin', 0, 'Zmax', 1);
stop_conditions = struct('maxFlips', test_params.max_flips, 'tolerance', test_params.tolerance);
start_conditions = struct('start', test_params.start_point, 'length', test_params.length);
obj = NelderMeadMethod(fobj, bounds, stop_conditions, start_conditions, settings, range);
% plot ideal minimum
plot3(X(1), X(2), X(3), 'x', 'color', 'g', 'lineWidth', 1.5);