classdef TestWavefunction < matlab.unittest.TestCase
    % TESTWAVEFUNCTION Unit tests for the Wavefunction class using a mock spatial grid

    properties
        spatialGrid  % Mock spatial grid for testing
        wavefunction % Instance of the Wavefunction class for testing
    end

    methods (TestMethodSetup)
        function setupMockGridAndWavefunction(testCase)
            % Setup mock spatial grid and Wavefunction instance before each test
            leftBoundary = -10;
            rightBoundary = 10;
            numPoints = 2000;
            L = rightBoundary - leftBoundary;

            % Generate the spatial grid (consistent with GridManager)
            testCase.spatialGrid = leftBoundary + L * (0:numPoints-1) / numPoints;

            % Inject mock grid into the Wavefunction class
            testCase.wavefunction = Wavefunction(struct('spatialGrid', testCase.spatialGrid));
        end
    end

    methods (Test)
        function testGaussianInitialization(testCase)
            % Test Gaussian wavepacket initialization
            center = 0; % Center of the Gaussian
            width = 1;  % Width of the Gaussian

            % Initialize the Gaussian wavepacket
            testCase.wavefunction = testCase.wavefunction.initializeGaussian(center, width);

            % Access the initialized state
            state = testCase.wavefunction.state;
            testCase.assertNotEmpty(state, 'Wavefunction state is empty.');

            % Check normalization (sum of probabilities should equal 1)
            normWave = sum(abs(state).^2);
            testCase.verifyEqual(normWave, 1, 'AbsTol', 1e-6, ...
                'Gaussian wavepacket is not normalized.');

            % Check peak position
            [~, peakIndex] = max(abs(state));
            peakPosition = testCase.spatialGrid(peakIndex);
            testCase.verifyEqual(peakPosition, center, 'AbsTol', 1e-6, ...
                'Gaussian wavepacket peak is not at the correct position.');
        end

        function testHermiteBasisGeneration(testCase)
            % Test Hermite basis generation
            n = 2; % Quantum number
            center = 0;
            M = 1; % Mass
            OMEGA = 1; % Angular frequency
            HBAR = 1; % Reduced Planck constant

            % Generate Hermite basis
            [basis, energy] = testCase.wavefunction.generateHermiteBasis(n, center, M, OMEGA, HBAR);

            % Check normalization (sum of probabilities should equal 1)
            normBasis = sum(abs(basis).^2);
            testCase.verifyEqual(normBasis, 1, 'AbsTol', 1e-6, ...
                'Hermite basis is not normalized.');

            % Check energy eigenvalue
            expectedEnergy = HBAR * OMEGA * (n + 0.5);
            testCase.verifyEqual(energy, expectedEnergy, 'AbsTol', 1e-6, ...
                'Energy eigenvalue of Hermite basis is incorrect.');
        end

        function testProbabilityDensityCalculation(testCase)
            % Test probability density calculation
            center = 0;
            width = 1;

            % Initialize the Gaussian wavepacket
            testCase.wavefunction = testCase.wavefunction.initializeGaussian(center, width);

            % Compute probability density
            probDensity = testCase.wavefunction.computeProbabilityDensity();
            testCase.assertNotEmpty(probDensity, 'Probability density is empty.');

            % Check normalization of probability density (should equal 1)
            integratedProb = sum(probDensity);
            testCase.verifyEqual(integratedProb, 1, 'AbsTol', 1e-6, ...
                'Probability density does not integrate to 1.');
        end

        function testStateNormalization(testCase)
            % Test the normalization function
            state = cos(pi * testCase.spatialGrid / 10); % Example arbitrary state
            normalizedState = testCase.wavefunction.normalize(state);

            % Check normalization (sum of probabilities should equal 1)
            normState = sum(abs(normalizedState).^2);
            testCase.verifyEqual(normState, 1, 'AbsTol', 1e-6, ...
                'State normalization is incorrect.');
        end
    end
end
