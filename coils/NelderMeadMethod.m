classdef NelderMeadMethod < handle

properties (Access = private)
    settings = struct('range', 4, 'step', 0.01, 'slices', 10, 'func_counter', 0, 'dimension', 3);
    stop_conditions = struct('maxFlips', 100, 'minHalving', 1e-2, 'minMargin', 1e-5);
    start_conditions = struct('start', [0 0 0], 'length', 0.25);
    f_objective = [];
    bounds = {};
    polytope = {};
    result = struct('halvings', 0, 'flips', 0);
end
    
methods
    % object constructor
    function this = NelderMeadMethod(f_objective, bounds, stop_conditions, start_conditions, settings)
        if settings.dimension ~= length(start_conditions.start)
            error("Input dimensions don't match");
        end
        
        this.f_objective = f_objective;
        this.bounds = bounds;
        this.stop_conditions = stop_conditions;
        this.start_conditions = start_conditions;
        this.settings = settings;
        
        % plot setup
        view(135,20)
        hold on
        colormap(hot)
        axis equal

        % setup first simplex
        this.polytope{end+1} = bsxfun(@plus, this.start_conditions.length.*this.simplexCoordinates(this.settings.dimension), this.start_conditions.start);
        
        % compute algorythm
        this.loop();
        
        % draw f objective
        this.plotFunction(this.f_objective, false);
        % draw bounds
        if ~isempty(this.bounds) > 0
            for k = 1:length(this.bounds)
                this.plotFunction(this.bounds{k}, true);
            end
        end
        % draw polytope
        for k = 1:length(this.polytope)
            this.drawSimplex(this.polytope{k}, [0.7 0.7 0.7])
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
        p = ones(1, length(V));
        if ~isempty(this.bounds)
            for j = 1:length(this.bounds)
                for k = 1:length(V)
                    if (this.bounds{j}(V(k, :)) < 0)
                        p(1, k) = p(1, k)*100;
                    end
                end
            end
        end
    end
    
    % returns the maximum vertex of the last simplex
    function i = findMaximumVertex(this, penality)
        fval = zeros(1, length(penality));
        for i = 1:length(penality)
            fval(i) = this.evaluate(this.f_objective, this.polytope{end}(i, :))*penality(i);
        end
        [~, i] = max(fval);
    end
    
    % returns the minimum vertex of the last simplex
    function i = findMinimumVertex(this, penality)
        fval = zeros(1, length(penality));
        for i = 1:length(penality)
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
        for i = 1:length(s)
           if i ~= j
              s(i, :) = (s(i, :)+s(j, :))/2; 
           end
        end
        this.result.halvings = this.result.halvings + 1;
    end
    
    function V = evaluate(this, f, v)
        V = f(v);
    end
    
    function V = clearBounds(this, V)
        % avoid cluttering the figure by plotting 0 values
        V(V == 0) = NaN;
        % zero the non-zero values to darken them in the heatmap
        V = V.*0;
    end
    
    function plotFunction(this, f, isBound)
        % divide plot in blocks
        if this.settings.dimension == 3
            X = -this.settings.range:this.settings.range/this.settings.slices:this.settings.range;
            Y = X;
            Z = X;
        elseif this.settings.dimension == 2
            [X, Y] = meshgrid(-this.settings.range:this.settings.step:this.settings.range);
        end
        len = length(X);
        
        % 3D plot
        if this.settings.dimension == 3
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
                    V = this.clearBounds(V);
                else
                    V
                    min(V)
                    max(V)
                end
                % properly offset the slice according to the number of function
                % passed via constructor
                K = ones(len, len)*(k-len/2)*(this.settings.range/this.settings.slices)+((this.settings.range/this.settings.slices)/(length(this.bounds)+length(this.f_objective)))*(this.settings.func_counter);
                surf(X, Y, K, 'CData', V, 'FaceAlpha', 0.3, 'LineStyle', 'none', 'FaceColor', 'interp');
            end
        elseif this.settings.dimension == 2
            V = zeros(len, len);
            for i = 1:len
                for j = 1:len
                    if isBound
                        V(i, j) = f([X(i) Y(j)]) < 0;
                    else
                        V(i, j) = f([X(i) Y(j)]);
                    end
                end
            end
            if isBound
                V = this.clearBounds(V);
            end
            surf(X, Y, V);
        end
        this.settings.func_counter = this.settings.func_counter + 1;
    end

    % computes the Cartesian coordinates of a simplex with centroid of 0
    function x = simplexCoordinates(this, n)
        x(1:n, 1:n+1) = 0.0;
        for j = 1:n
            x(j, j) = 1.0;
        end
        x(1:n, n+1) = (1.0-sqrt(1.0+n))/n;
        c(1:n, 1) = sum(x(1:n, 1:n+1), 2)/(n+1);
        for j = 1:n+1
            x(1:n, j) = x(1:n, j) - c(1:n, 1);
        end
        s = norm(x(1:n, 1));
        x(1:n, 1:n+1) = x(1:n, 1:n+1)/s;
        if this.settings.dimension == 3
            x = x./1.633;
        end
        x = x';
    end
    
    % draws a simplex
    function drawSimplex(this, P, color)
        if this.settings.dimension == 3
            f = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
        elseif this.settings.dimension == 2
            f = [1 2 3];
        end
        patch('Faces', f, 'Vertices', P, 'EdgeColor', 'w', 'FaceColor', color, 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    end
    
    % returns the centroid of a regular tetrahedron
    function x = centroid(this, t)
        if this.settings.dimension == 3
            x = [sum(t(:, 1))/length(t) sum(t(:, 2))/length(t) sum(t(:, 3))/length(t)];
        elseif this.settings.dimension == 2
            x = [sum(t(:, 1))/length(t) sum(t(:, 2))/length(t)];
        end
    end
    
    % returns the height of a regular tetrahedron
    function h = simplexHeight(this, P)
        if this.settings.dimension == 3
            h = sqrt(6)/3.*norm(P(1, :) - P(2, :));
        elseif this.settings.dimension == 2
            h = sqrt(3)/2.*norm(P(1, :) - P(2, :));
        end
    end
end
    
end