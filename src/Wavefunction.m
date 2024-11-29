classdef Wavefunction
    % WAVEFUNCTION Handles initialization, normalization, and utility functions.

    properties
        state       % Wavefunction state (complex array)
        grid        % Spatial grid associated with the wavefunction
        params      % Parameters for the wavefunction
    end

    methods
        %% Constructor
        function obj = Wavefunction(gridManager)
            % Constructor to initialize with spatial grid
            obj.grid = gridManager.getSpatialGrid();
        end

        %% Initialization
        % function obj = initializeGaussian(obj, center, width)
        %     % Initialize a Gaussian wavepacket
        %     psi = exp(-(obj.grid - center).^2 / (2 * width^2));
        %     obj.state = obj.normalize(psi);
        % end

        function obj = initializeHarmonicOscillator(obj, params)
            % Initialize superposition of SHO states
            obj.params = params;
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
            % normalization = (M * OMEGA / pi / HBAR)^(1/4) / sqrt(2^n * factorial(n));  % Normalization constant
            basis = 1 * H_n .* exp(-M * OMEGA* x.^2 / (2 * HBAR));
            basis = obj.normalize(basis);
            energy = HBAR * OMEGA * (n + 0.5);
        end

        %% Inner Product
        function innerProd = innerProduct(obj, otherKet)
            % INNERPRODUCT Computes the inner product with another ket
            % INPUT:
            %   - otherKet: The wavefunction (ket) to take the inner product with
            % OUTPUT:
            %   - innerProd: The resulting complex scalar value of the inner product
    
            % Validate dimensions
            if length(obj.state) ~= length(otherKet)
                error('Wavefunction:DimensionMismatch', ...
                      'The input ket must have the same dimension as the current wavefunction.');
            end
    
            % Compute the inner product using the DVR basis
            innerProd = sum(conj(obj.state) .* otherKet);
        end

        %% Utility Functions
        function normalizedState = normalize(obj, state)
            % Normalize a wavefunction
            normalizedState = state / sqrt(sum(abs(state).^2));
            %disp(sum(abs(normalizedState).^2));
        end

        function probDensity = computeProbabilityDensity(obj)
            % Compute the probability density of the wavefunction
            probDensity = abs(obj.state).^2;
        end
    end
end
