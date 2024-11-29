classdef TestGridManager < matlab.unittest.TestCase
    % TESTGRIDMANAGER Unit tests for the GridManager class

    methods (Test)
        function testGridGeneration(testCase)
            % Test generation of spatial and momentum grids
            leftBoundary = -10;
            rightBoundary = 10;
            numPoints = 2000;

            % Instantiate the GridManager
            gridManager = GridManager(leftBoundary, rightBoundary, numPoints);

            % Verify spatial grid properties
            spatialGrid = gridManager.getSpatialGrid();
        end

        function testMomentumGridSymmetry(testCase)
            % Test symmetry of the momentum grid
            leftBoundary = -10;
            rightBoundary = 10;
            numPoints = 2000;

            % Instantiate the GridManager
            gridManager = GridManager(leftBoundary, rightBoundary, numPoints);

            % Verify momentum grid size
            momentumGrid = gridManager.getMomentumGrid();
            testCase.verifyEqual(length(momentumGrid), numPoints, ...
                'Momentum grid does not have the correct number of points.');

            % Verify momentum grid difference
            diffMomentum = abs(momentumGrid(2) - momentumGrid(1));
            testCase.verifyEqual(diffMomentum, 2 * pi / (rightBoundary-leftBoundary), 'AbsTol', 1e-6, ...
                'Momentum grid is interval is not 2*pi/L.');
        end

        function testGridBoundaries(testCase)
            % Test spatial grid boundaries
            leftBoundary = 5;
            rightBoundary = 15;
            numPoints = 100;

            % Instantiate the GridManager
            gridManager = GridManager(leftBoundary, rightBoundary, numPoints);

            % Verify spatial grid boundaries
            spatialGrid = gridManager.getSpatialGrid();
            testCase.verifyEqual(spatialGrid(1), leftBoundary, 'AbsTol', 1e-6, ...
                'Spatial grid does not start at the left boundary.');
            
            % Verify size of spatial grid
            testCase.verifyEqual(length(spatialGrid), numPoints, ...
                'Spatial grid does not have the correct number of points.');
        end
    end
end