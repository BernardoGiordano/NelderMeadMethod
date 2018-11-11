classdef NelderMeadMethod < handle

properties (Access = private)
    f_objective = [];
    bounds = {};
    stop_conditions = struct('steps', 100, 'minArea', 1e-5, 'minHalving', 1e-2, 'minMargin', 1e-5);
    start_conditions = struct('start', [0 0 0], 'area', 1);
    range = 1;
    slices = 20;
    func_counter = 0;
    polytope = {};
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
        colormap(hot)
        axis equal

        % setup first simplex
        l = sqrt(4*this.start_conditions.area/sqrt(3));
        this.polytope{end+1} = bsxfun(@plus, l.*this.simplexCoordinates(3)./1.633, this.start_conditions.start);
        
        % draw f objective
        this.plotFunction(this.f_objective, false);
        % draw bounds
        for k = 1:length(this.bounds)
            this.plotFunction(this.bounds{k}, true);
        end
        % draw polytope
        for k = 1:length(this.polytope)
            this.drawSimplex(this.polytope{k}, 'white')
        end
    end
    
    function V = evaluate(this, f, v)
        V = f(v);
    end
    
    function plotFunction(this, f, isBound)
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

    % computes the Cartesian coordinates of a simplex with centroid of 0
    function x = simplexCoordinates(this, n)
        x(1:n, 1:n+1) = 0.0;
        for j = 1:n
            x(j, j) = 1.0;
        end
        a = (1.0-sqrt(1.0+n))/n;
        x(1:n, n+1) = a;
        c(1:n, 1) = sum(x(1:n, 1:n+1), 2)/(n+1);
        for j = 1:n+1
            x(1:n, j) = x(1:n, j) - c(1:n, 1);
        end
        s = norm(x(1:n, 1));
        x(1:n, 1:n+1) = x(1:n, 1:n+1)/s;
        x = x';
    end
    
    % draws a simplex
    function drawSimplex(this, P, color)
        f = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
        patch('Faces', f ,'Vertices', P(:, 1:3), 'EdgeColor', 'k', 'FaceColor', color, 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    end
end
    
end