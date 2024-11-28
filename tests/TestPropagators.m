classdef TestPropagators < matlab.unittest.TestCase
    % TestPropagators Unit tests for the Propagators class.
    
    properties
        propagator % Instance of the Propagators class
        gridManager % Example grid manager
        spatialGrid % Spatial grid
        momentumGrid % Momentum grid
        dt % Time step
        totalTime % Total simulation time
        numSteps % Number of time steps
    end
    
    methods (TestMethodSetup)
        function createTestInstance(testCase)
            % Define test parameters
            testCase.gridManager = GridManager(100, 100, 10000); % Example grid manager
            testCase.spatialGrid = testCase.gridManager.spatialGrid;
            testCase.momentumGrid = testCase.gridManager.momentumGrid;
            testCase.totalTime = 2 * pi; % Total time
            testCase.numSteps = 20000; % Number of steps
            testCase.dt = 2 * pi / 20000; % Time step
            
            % Define test energy functions
            kineticEnergyFunc = @(p) p.^2 / 2; % Kinetic energy (free particle)
            potentialEnergyFunc = @(x, t) x.^2 + sin(t); % Harmonic potential with time-dependence
            
            % Create instance of the Propagators class
            testCase.propagator = Propagators(kineticEnergyFunc, potentialEnergyFunc);
        end
    end
    
    methods (Test)
        function testKineticPropagator(testCase)
            % Test the kinetic propagator generation
            UT = testCase.propagator.generateKineticPropagator(testCase.momentumGrid, testCase.dt);
            
            % Verify size and type
            testCase.verifySize(UT, size(testCase.momentumGrid));
            testCase.verifyClass(UT, 'double');
            
            % Verify correctness for a simple case
            p = 1; % Example momentum
            expected = exp(-1i * (p^2 / 2) * testCase.dt);
            actual = exp(-1i * testCase.propagator.kineticEnergyFunc(p) * testCase.dt);
            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-12);
        end
        
        function testPotentialPropagators(testCase)
            % Test the potential propagators generation
            UV_list = testCase.propagator.generatePotentialPropagators( ...
                testCase.spatialGrid, testCase.dt, testCase.totalTime, testCase.numSteps);
            
            % Verify size and type of the output
            testCase.verifyEqual(numel(UV_list), testCase.numSteps);
            testCase.verifyClass(UV_list, 'cell');
            
            % Verify correctness for the first time step
            t = 0; % First time step
            potentialEnergy = testCase.propagator.potentialEnergyFunc(testCase.spatialGrid, t);
            expected = exp(-1i * potentialEnergy * testCase.dt / 2);
            testCase.verifyEqual(UV_list{1}, expected, 'AbsTol', 1e-12);
        end
    end
end
