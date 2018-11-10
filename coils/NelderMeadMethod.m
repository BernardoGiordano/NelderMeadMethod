classdef NelderMeadMethod < handle

properties (Access = private)
    f_objective = [];
    bounds = {};
    stop_conditions = struct('steps', 100, 'minArea', 1e-5, 'minHalving', 1e-2, 'minMargin', 1e-5);
    start_conditions = struct('start', [0,0,0], 'area', 1);
    range = 1;
    slices = 20;
    func_counter = 0;
end
    
methods
    % object constructor
    function this = NelderMeadMethod(f_objective, bounds, stop_conditions, start_conditions)
        this.f_objective = f_objective;
        this.bounds = bounds;
        this.stop_conditions = stop_conditions;
        this.start_conditions = start_conditions;
        
        % plot setup
        view(3)
        hold on
        colormap(hot);
        
        % draw f objective
        this.plot_function(this.f_objective, false);
        % draw bounds
        for k = 1:length(this.bounds)
            this.plot_function(this.bounds{k}, true);
        end
    end
    
    function plot_function(this, f, isBound)
        % divide plot in blocks
        X = -this.range:this.range/this.slices:this.range;
        Y = X;
        Z = X;
        len = length(X);
        
        for k = 1:len
            V = zeros(len, len);
            for i = 1:len
                for j = 1:len
                    if isBound
                        V(i, j) = f([X(i) Y(j) Z(k)]) < 0;
                    else
                        V(i, j) = f([X(i) Y(j) Z(k)]);
                    end
                end
            end
            if isBound
                % avoid cluttering the figure by plotting 0 values
                V(V == 0) = NaN;
                % zero the non-zero values to darken them in the heatmap
                V = V.*0;
            end
            % properly offset the slice according to the number of function
            % passed via constructor
            K = ones(len, len)*(k-len/2)*(this.range/this.slices)+((this.range/this.slices)/(length(this.bounds)+length(this.f_objective)))*(this.func_counter);
            surf(X, Y, K, 'CData', V, 'FaceAlpha', 0.3, 'LineStyle', 'none', 'FaceColor', 'interp');
        end
        this.func_counter = this.func_counter + 1;
    end
end
    
end