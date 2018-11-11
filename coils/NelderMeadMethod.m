classdef NelderMeadMethod < handle

properties (Access = private)
    f_objective = [];
    bounds = {};
    stop_conditions = struct('maxFlips', 100, 'minArea', 1e-5, 'minHalving', 1e-2, 'minMargin', 1e-5);
    start_conditions = struct('start', [0 0 0], 'area', 1);
    range = 1;
    slices = 20;
    func_counter = 0;
    polytope = {};
    result = struct('halvings', 0, 'flips', 0)
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
        this.polytope{end+1} = bsxfun(@plus, l.*this.simplexCoordinates(3), this.start_conditions.start);
        
        % compute algorythm
        this.loop();
        
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
        % display results
        disp(this.result);
    end
    
    function loop(this)
        % 1) set penality for vertices out of bounds
        % 2) test flipping loop condition
        % 2.yes.1) find minimum
        % 2.yes.2) halve
        % 2.no.1) find maximum
        % 2.no.2 flip
        % 3) detect stop condition
        shouldContinue = true;
        while shouldContinue
            penality = this.getPenality(this.polytope{end});
            if this.isFlippingForever()
                j = this.findMinimumVertex(penality);
                S = this.halve(j);
                this.polytope{end+1} = S;
            else
                j = this.findMaximumVertex(penality);
                S = this.flip(j);
                this.polytope{end+1} = S;
            end
            
            % check stop conditions
            if this.result.flips >= this.stop_conditions.maxFlips
                % TODO: check more stop conditions
                shouldContinue = false;
            end
        end
    end
    
    % returns if flipping loop condition occurred
    function v = isFlippingForever(this)
        v = false;
        n = 8;
        if length(this.polytope) > n
            subset = this.polytope(end-n+1:end);
            for j = 1:n
                for k = 1:n
                    if j ~= k
                        if cell2mat(subset(k)) == cell2mat(subset(j))
                            v = true;
                            return
                        end
                    end
                end
            end
        end
    end
    
    % returns the penality array for a given simplex
    function p = getPenality(this, V)
        p = ones(1, 4);
        for j = length(this.bounds)
            for k = length(V)
                if (this.bounds{j}(V(k, :)) < 0)
                    p(1, k) = p(1, k)*100;
                end
            end
        end
    end
    
    % returns the maximum vertex of the last simplex
    function i = findMaximumVertex(this, penality)
        fval = zeros(1, 4);
        for i = 1:4
            fval(i) = this.evaluate(this.f_objective, this.polytope{end}(i, :))*penality(i);
        end
        [~, i] = max(fval);
    end
    
    % returns the minimum vertex of the last simplex
    function i = findMinimumVertex(this, penality)
        fval = zeros(1, 4);
        for i = 1:4
            fval(i) = this.evaluate(this.f_objective, this.polytope{end}(i, :))*penality(i);
        end
        [~, i] = min(fval);
    end
    
    % flips the j-th vertex of the last simplex in the polytope array 
    function s = flip(this, j)
        s = this.polytope{end};
        c = this.centroid(s);
        h = this.simplexHeight(s);
        % find direction of jth_vertex - centroid vector
        v = c - s(j, :);
        % get versor of said vector
        v = v/norm(v);
        % jth' vertex = jth vertex + 2*h*versor
        s(j, :) = s(j, :) + 2*h.*v;
        % increment flips counter
        this.result.flips = this.result.flips + 1;
    end
    
    % halve the last simplex keeping the j-th vertex
    function s = halve(this, j)
        s = this.polytope{end};
        for i = 1:4
           if i ~= j
              s(i, :) = (s(i, :)+s(j, :))/2; 
           end
        end
        this.result.halvings = this.result.halvings + 1;
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
        x = x'./1.633;
    end
    
    % draws a simplex
    function drawSimplex(this, P, color)
        f = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
        patch('Faces', f ,'Vertices', P(:, 1:3), 'EdgeColor', 'k', 'FaceColor', color, 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    end
    
    % returns the centroid of a regular tetrahedron
    function x = centroid(this, t)
        x = [sum(t(:, 1))/4 sum(t(:, 2))/4 sum(t(:, 3))/4];
    end
    
    % returns the height of a regular tetrahedron
    function h = simplexHeight(this, P)
        h = sqrt(6)/3.*norm(P(1, :) - P(2, :));
    end
end
    
end