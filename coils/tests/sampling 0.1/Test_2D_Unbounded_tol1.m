addpath(genpath('../'));

clear all
close all
clc

global hZ
% hidden parameters
hI = 5;
hR = 0.8;
hZ = 0.7;
X = [hI hR];

disp("Test parameters")
%%%%%%%%%%%%%% EDIT ME %%%%%%%%%%%%%%
test_params = struct( ...
    'dimension', 2, ...
    'plot', true, ...
    'minimum', X, ...
    'sampling_step', 0.1, ...
    'sampling_range', 2, ...
    'sample_amount', -1, ...
    'max_flips', 1000, ...
    'tolerance', 0.1, ...
    'minLength', 1e-3, ...
    'start_point', [3 0.3], ...
    'length', 0.3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = -test_params.sampling_range:test_params.sampling_step:test_params.sampling_range;
test_params.sample_amount = length(z);
disp(test_params)

%%%%%%%%%%%%%% EDIT ME %%%%%%%%%%%%%%
bounds = {};
fobj = @(x)((1/mean(Bz2d(X, z)))*(norm(Bz2d(x, z) - Bz2d(X, z))*sqrt(inv(length(z)))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simplex algorythm
settings = struct('step', test_params.sampling_step, 'slices', length(z), 'plot', test_params.plot, 'dimension', test_params.dimension);
range = struct('Xmin', 0, 'Xmax', 6, 'Ymin', 0, 'Ymax', 1);
stop_conditions = struct('maxFlips', test_params.max_flips, 'tolerance', test_params.tolerance, 'minLength', test_params.minLength);
start_conditions = struct('start', test_params.start_point, 'length', test_params.length);
obj = NelderMeadMethod(fobj, bounds, stop_conditions, start_conditions, settings, range);
if test_params.plot
    % plot ideal minimum
    plot(X(1), X(2), 'x', 'color', 'g', 'lineWidth', 1.5);
end

function res = Bz2d(X, z)
    global hZ;
    I = X(1);
    R = X(2);
    Z = hZ;
    res = Bz([I R Z], z);
end