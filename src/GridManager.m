classdef GridManager
    % GRIDMANAGER Handles spatial and momentum grid generation.

    properties
        leftBoundary   % Start of spatial domain
        rightBoundary  % End of spatial domain
        numPoints      % Number of spatial grid points
        spatialGrid    % Spatial grid array
        momentumGrid   % Momentum grid array
    end

    methods
        function obj = GridManager(leftBoundary, rightBoundary, numPoints)
            % Constructor to initialize grid parameters
            obj.leftBoundary = leftBoundary;
            obj.rightBoundary = rightBoundary;
            obj.numPoints = numPoints;
            [obj.spatialGrid, obj.momentumGrid] = obj.generateGrids();
        end
    
        function [spatialGrid, momentumGrid] = generateGrids(obj)
            % Generate spatial and momentum grids
            L = obj.rightBoundary - obj.leftBoundary;
            spatialGrid = obj.leftBoundary + L * (0:obj.numPoints-1) / obj.numPoints;
            momentumGrid = (2 * pi / L) * [0:floor(obj.numPoints/2)-1, -ceil(obj.numPoints/2):-1];
        end
    end
end
