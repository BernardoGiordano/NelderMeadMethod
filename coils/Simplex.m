classdef Simplex

properties
    f_objective = 0;
    bounds = [];
    stop_conditions = struct ('step', 100, 'minArea', 1e-5, 'minHalving', 1e-2, 'minMargin', 1e-5);
    start_conditions = struct ('start', [0,0,0], 'area', 1);
    % si definisce prima il campo e poi il contenuto, es. ('step',[])
    % [] = contenuto vuoto
    % stop_conditions(1).step = 100; attribuisce un valore nella posizione 1 
    % stop_conditions(1).step richiama l'elemento 1  
    field = 1;
    slices = 10;
end
    
methods
    % object constructor
    function obj = Simplex(f_objective, bounds, stop_conditions, start_conditions)
        obj.f_objective = f_objective;
        obj.bounds = bounds;
        obj.stop_conditions = stop_conditions;
        obj.start_conditions = start_conditions;
        % plot setup
        view(3)
        hold on
        % draw f objective
        % obj.draw(f_objective);
        % draw bounds
        for bound = obj.bounds
            obj.draw(bound);
        end
    end
    
    %draws a generic 3D function as a heatmap 
    function draw(obj, f)
        X = -obj.field:obj.field/obj.slices:obj.field;
        Y = -obj.field:obj.field/obj.slices:obj.field;
        Z = -obj.field:obj.field/obj.slices:obj.field;
        V = zeros(obj.slices, length(X), length(Y));

        %calculates the function values in a 3-dimensional grid
        for k = 1:length(Z)
            for i = 1:length(Y)
                for j = 1:length(X)
                    V(k, i, j) = f([X(j) Y(i) Z(k)]);
                end
            end
        end

        for j = 1:length(Z)
            %calculate the plane
            offset = ((j-(length(Z)/2))*(obj.field/obj.slices));
            Z = ones(length(X), length(Y))*offset;

            %draw the plane with the iso-level colors of the function
            %computed in its coordinates
            colormap(hot);
            surf(X, Y, Z, 'CData', reshape(V(j,:,:), [length(X), length(Y)]), 'FaceAlpha', 0.3, 'LineStyle', 'none', 'FaceColor', 'interp');
        end
    end
end
    
end