classdef NelderMeadMethod < handle

properties (Access = private)
    settings = struct('range', 4, 'step', 0.01, 'slices', 10, 'dimension', 3);
    range = struct('Xmin', 0, 'Xmax', 6, 'Ymin', 0, 'Ymax', 1, 'Zmin', 0, 'Zmax', 1);
    stop_conditions = struct('maxFlips', 100, 'minHalving', 1e-2, 'minMargin', 1e-5);
    start_conditions = struct('start', [0 0 0], 'length', 0.25);
    internal = struct('func_counter', 0);
    f_objective = [];
    bounds = {};
    polytope = {};
    result = struct('halvings', 0, 'flips', 0);
end
    
methods
    % object constructor
    function this = NelderMeadMethod(f_objective, bounds, stop_conditions, start_conditions, settings, range)
        if settings.dimension ~= length(start_conditions.start)
            error("Input dimensions don't match");
        end
        
        this.f_objective = f_objective;
        this.bounds = bounds;
        this.stop_conditions = stop_conditions;
        this.start_conditions = start_conditions;
        this.settings = settings;
        this.range = range;
        
        % plot setup
        if this.settings.dimension == 3
            view(135,20)
            axis([this.range.Xmin this.range.Xmax this.range.Ymin this.range.Ymax this.range.Zmin this.range.Zmax])
        elseif this.settings.dimension == 2
            % TODO
            axis([this.range.Xmin this.range.Xmax this.range.Ymin this.range.Ymax])
        end
        hold on
        colormap(hot)
        axis equal

        % setup first simplex
        this.polytope{end+1} = bsxfun(@plus, this.start_conditions.length.*this.simplexCoordinates(this.settings.dimension), [this.start_conditions.start, 0]);

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
        shouldContinue = true;
        while shouldContinue
            s = this.polytope{end};
            % set penality for out of bounds vertices
            penality = this.getPenality(s);
            % age counter is too big for at least one of the vertices
            if this.isTooOld(s)
                % halve with the minimum vertex as pivot
                j = this.findMinimumVertex(s, penality);
                new_s = this.halve(s, j);
            else
                % flip with respect to the maximum vertex
                j = this.findMaximumVertex(s, penality);
                new_s = this.flip(j);
                % if the new maximum value is bigger than the previous one after
                % flipping, flip with respect of the 2th maximum vertex
                v_s = this.evaluate(this.f_objective, s(j, 1:this.settings.dimension));
                v_new_s = this.evaluate(this.f_objective, new_s(j, 1:this.settings.dimension));
                if v_s < v_new_s
                    new_s = this.flip(this.find2thMaximumVertex(s, penality));
                end
                % check for out of range conditions
                if this.isOutOfRange(new_s)
                    new_s = this.halve(s, this.findMinimumVertex(s, penality));
                end
            end
            % append next simplex
            this.polytope{end+1} = new_s;
            % check stop conditions
            if this.result.flips >= this.stop_conditions.maxFlips
                % TODO: check more stop conditions
                shouldContinue = false;
            end
        end
    end
    
    % true if a simplex vertex is out of range
    function r = isVertexOutOfRange(this, s, i)
        r = false;
        [~, M] = size(s);
        for j = 1:M-1
            if s(i, j) > this.settings.range || s(i, j) < -this.settings.range
                r = true;
                break
            end
        end
    end
    
    % true if a simplex is out of range
    function r = isOutOfRange(this, s)
        r = false;
        [N, ~] = size(s);
        for i = 1:N
            if this.isVertexOutOfRange(s, i)
                r = true;
                break
            end
        end
    end
    
    % returns if the age counter is more than 
    function v = isTooOld(this, s)
        v = false;
        % TODO: parametrize n
        n = 7;
        [N, ~] = size(s);
        for i = 1:N
            if s(i, this.settings.dimension + 1) == n
                v = true;
                break
            end
        end
    end
    
    % returns the penality array for a given simplex
    function p = getPenality(this, V)
        p = ones(1, length(V));
        if ~isempty(this.bounds)
            for j = 1:length(this.bounds)
                for k = 1:length(V)
                    if (this.bounds{j}(V(k, 1:this.settings.dimension)) < 0)
                        p(1, k) = p(1, k)*100;
                    end
                end
            end
        end
    end
    
    % returns the maximum vertex of the last simplex
    function i = findMaximumVertex(this, s, penality)
        fval = zeros(1, length(penality));
        for i = 1:length(penality)
            if this.isVertexOutOfRange(s, i)
                fval(i) = NaN;
            else
                fval(i) = this.evaluate(this.f_objective, s(i, 1:this.settings.dimension))*penality(i);
            end
        end
        [~, i] = max(fval);
    end
    
    % returns the minimum vertex of the last simplex
    function i = findMinimumVertex(this, s, penality)
        fval = zeros(1, length(penality));
        for i = 1:length(penality)
            if this.isVertexOutOfRange(s, i)
                fval(i) = NaN;
            else
                fval(i) = this.evaluate(this.f_objective, s(i, 1:this.settings.dimension))*penality(i);
            end
        end
        [~, i] = min(fval);
    end
    
    % returns the 2th maximum vertex of the last simplex
    function i = find2thMaximumVertex(this, s, penality)
        fval = zeros(1, length(penality));
        for i = 1:length(penality)
            fval(i) = this.evaluate(this.f_objective, s(i, 1:this.settings.dimension))*penality(i);
        end
        [~, i] = max(fval);
        fval(1, i) = NaN;
        [~, i] = max(fval);
    end
    
    % flips the j-th vertex of the last simplex in the polytope array 
    function s = flip(this, j)
        s = this.polytope{end};
        c = this.centroid(s);
        h = this.simplexHeight(s);
        % find direction of jth_vertex - centroid vector
        v = c - s(j, 1:this.settings.dimension);
        % get versor of said vector
        v = v/norm(v);
        % jth' vertex = jth vertex + 2*h*versor
        s(j, 1:this.settings.dimension) = s(j, 1:this.settings.dimension) + 2*h.*v;
        % increment flips counter
        this.result.flips = this.result.flips + 1;
        % increment age counter
        s(:, this.settings.dimension + 1) = s(:, this.settings.dimension + 1) + 1;
        s(j, this.settings.dimension + 1) = 0;
    end
    
    % halve the last simplex keeping the j-th vertex
    function s = halve(this, s, j)
        for i = 1:length(s)
           if i ~= j
              s(i, 1:this.settings.dimension) = (s(i, 1:this.settings.dimension) + s(j, 1:this.settings.dimension))/2; 
           end
        end
        % increment halvings counter
        this.result.halvings = this.result.halvings + 1;
        % set age counter to 0 for each vertex
        s(:, this.settings.dimension + 1) = 0;
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
                end
                % properly offset the slice according to the number of function
                % passed via constructor
                K = ones(len, len)*(k-len/2)*(this.settings.range/this.settings.slices)+((this.settings.range/this.settings.slices)/(length(this.bounds)+length(this.f_objective)))*(this.internal.func_counter);
                surf(X, Y, K, 'CData', V, 'FaceAlpha', 0.3, 'LineStyle', 'none', 'FaceColor', 'interp');
            end
        % 2D plot
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
        this.internal.func_counter = this.internal.func_counter + 1;
    end

    % computes the Cartesian coordinates of a simplex with centroid in 0
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
        x = cat(2, x', zeros(length(x), 1));
    end
    
    % draws a simplex
    function drawSimplex(this, P, color)
        if this.settings.dimension == 3
            f = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
        elseif this.settings.dimension == 2
            f = [1 2 3];
        end
        patch('Faces', f, 'Vertices', P(:, 1:this.settings.dimension), 'EdgeColor', 'w', 'FaceColor', color, 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    end
    
    % returns the centroid of a simplex
    function x = centroid(this, t)
        if this.settings.dimension == 3
            x = [sum(t(:, 1))/length(t) sum(t(:, 2))/length(t) sum(t(:, 3))/length(t)];
        elseif this.settings.dimension == 2
            x = [sum(t(:, 1))/length(t) sum(t(:, 2))/length(t)];
        end
    end
    
    % returns the height of a simplex
    function h = simplexHeight(this, P)
        if this.settings.dimension == 3
            h = sqrt(6)/3.*norm(P(1, 1:this.settings.dimension) - P(2, 1:this.settings.dimension));
        elseif this.settings.dimension == 2
            h = sqrt(3)/2.*norm(P(1, 1:this.settings.dimension) - P(2, 1:this.settings.dimension));
        end
    end
end
    
end