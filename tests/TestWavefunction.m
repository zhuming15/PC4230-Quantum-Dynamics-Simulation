classdef TestWavefunction < matlab.unittest.TestCase
    % TESTWAVEFUNCTION Unit tests for the Wavefunction class using DVR basis
    
    properties
        spatialGrid  % Spatial grid for testing
        wavefunction % Instance of the Wavefunction class for testing
    end

    methods (TestMethodSetup)
        function setupGridAndWavefunction(testCase)
            % Setup spatial grid and Wavefunction instance before each test
            testCase.spatialGrid = linspace(-10, 10, 2000); % Example grid
            testCase.wavefunction = Wavefunction(testCase.spatialGrid); % Instantiate Wavefunction
        end
    end

    methods (Test)
        function testGaussianInitialization(testCase)
            % Test Gaussian wavepacket initialization
            center = 0; % Center of the Gaussian
            width = 1;  % Width of the Gaussian

            % Explicitly reassign the updated object
            testCase.wavefunction = testCase.wavefunction.initializeGaussian(center, width);

            % Access the initialized state
            state = testCase.wavefunction.state;
            testCase.assertNotEmpty(state, 'Wavefunction state is empty.');

            % Check normalization in DVR (sum should equal 1)
            normWave = sum(abs(state).^2);
            testCase.verifyEqual(normWave, 1, 'AbsTol', 1e-6, ...
                'Gaussian wavepacket is not normalized.');
        end

        function testProbabilityDensityCalculation(testCase)
            % Test probability density calculation
            center = 0;
            width = 1;

            % Explicitly reassign the updated object
            testCase.wavefunction = testCase.wavefunction.initializeGaussian(center, width);

            % Compute probability density
            probDensity = testCase.wavefunction.computeProbabilityDensity();
            testCase.assertNotEmpty(probDensity, 'Probability density is empty.');

            % Check normalization of probability density in DVR
            integratedProb = sum(probDensity); % Should equal 1
            testCase.verifyEqual(integratedProb, 1, 'AbsTol', 1e-6, ...
                'Probability density does not integrate to 1.');
        end

        function testStateNormalization(testCase)
            % Test the normalization function
            state = cos(pi * testCase.spatialGrid / 10); % Example arbitrary state
            normalizedState = testCase.wavefunction.normalize(state);

            % Check normalization in DVR
            normState = sum(abs(normalizedState).^2);
            testCase.verifyEqual(normState, 1, 'AbsTol', 1e-6, ...
                'State normalization is incorrect.');
        end
    end
end
