classdef Wavefunction
    % WAVEFUNCTION Handles initialization, normalization, and utility functions.

    properties
        state       % Wavefunction state (complex array)
        grid        % Spatial grid associated with the wavefunction
    end

    methods
        %% Constructor
        function obj = Wavefunction(grid)
            % Constructor to initialize with spatial grid
            obj.grid = grid;
        end

        %% Initialization
        function obj = initializeGaussian(obj, center, width)
            % Initialize a Gaussian wavepacket
            psi = exp(-(obj.grid - center).^2 / (2 * width^2));
            obj.state = obj.normalize(psi);
        end

        function obj = initializeHarmonicOscillator(obj, params)
            % Initialize superposition of SHO states
            combinedState = zeros(size(obj.grid));
            for k = 1:numel(params.bases)
                [basis, ~] = obj.generateHermiteBasis(params.bases(k), params.center, ...
                    params.M, params.OMEGA, params.HBAR);
                combinedState = combinedState + params.coefficients(k) * basis;
            end
            obj.state = obj.normalize(combinedState);
        end

        %% Hermite Basis
        function [basis, energy] = generateHermiteBasis(obj, n, center, M, OMEGA, HBAR)
            x = obj.grid - center;
            xi = sqrt(M * OMEGA / HBAR) * x;
            H_n = hermiteH(n, xi);
            basis = H_n .* exp(-xi.^2 / 2);
            basis = obj.normalize(basis);
            energy = HBAR * OMEGA * (n + 0.5);
        end

        %% Utility Functions
        function normalizedState = normalize(~, state)
            % Normalize a wavefunction
            normalizedState = state / norm(state);
        end

        function probDensity = computeProbabilityDensity(obj)
            % Compute the probability density of the wavefunction
            probDensity = abs(obj.state).^2;
        end
    end
end
