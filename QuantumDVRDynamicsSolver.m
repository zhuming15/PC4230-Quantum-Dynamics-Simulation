classdef QuantumDVRDynamicsSolver
    properties
        leftBoundary    % Start of the spatial domain
        rightBoundary   % End of the spatial domain
        numPoints       % Number of grid points in the spatial domain
        spatialGrid     % Array of spatial grid points
        potentialFunc   % Handle to potential energy function V(x)
        hamiltonian     % Hamiltonian matrix in DVR
        eigenvalues     % Eigenvalues of the Hamiltonian
        eigenfunctions  % Eigenfunctions of the Hamiltonian
    end
    
    methods
        function obj = QuantumDVRDynamicsSolver(left, right, numPts)
            % Constructor to initialize the quantum dynamics solver
            obj.leftBoundary = left;
            obj.rightBoundary = right;
            obj.numPoints = numPts;
            obj.spatialGrid = linspace(left, right, numPts);
            obj.hamiltonian = [];  % Initialize to empty; will be calculated later
            obj.eigenvalues = [];
            obj.eigenfunctions = [];
        end
        
        function obj = generateHamiltonian(obj, potentialFunc)
            % Generates the Hamiltonian matrix using the DVR approach and specified potential
            h = (obj.rightBoundary - obj.leftBoundary) / (obj.numPoints - 1);  % Grid spacing
            obj.hamiltonian = zeros(obj.numPoints);
            kinetic = diag((pi^2 / (6 * h^2)) * (obj.numPoints^2 - 1) * ones(obj.numPoints, 1));  % Kinetic energy part
            potential = diag(arrayfun(potentialFunc, obj.spatialGrid));  % Potential energy part
            
            % Assembling the Hamiltonian matrix
            for i = 2:obj.numPoints
                for j = 1:i-1
                    kinetic(i, j) = (-1)^(i-j) * pi^2 / (3 * h^2 * sin(pi * (i-j) / obj.numPoints)^2);
                    kinetic(j, i) = kinetic(i, j);  % Hermitian matrix
                end
            end
            
            obj.hamiltonian = kinetic + potential;
        end
        
        function [obj, vec, val] = solveEigenproblems(obj)
            % Solves for eigenvalues and eigenfunctions of the Hamiltonian
            if isempty(obj.hamiltonian)
                error('Hamiltonian has not been generated. Call generateHamiltonian first.');
            end
            
            [vec, val] = eig(obj.hamiltonian);
            obj.eigenvalues = diag(val);
            obj.eigenfunctions = vec;
        end
        
        function [obj, figHandle] = plotEigenfunctions(obj, EigenfunctionsToPlot)
            % Plot a specified number of eigenfunctions
            if isempty(obj.eigenfunctions)
                error('Eigenfunctions have not been calculated. Call solveEigenproblems first.');
            end

            figHandle = figure;
            for k = EigenfunctionsToPlot
                plot(obj.spatialGrid, obj.eigenfunctions(:, k));
                hold on;
            end
            title('First few eigenfunctions of the quantum system');
            xlabel('Spatial coordinate');
            ylabel('Wavefunction amplitude');
            legend(arrayfun(@(x) ['\psi_', num2str(x)], EigenfunctionsToPlot, 'UniformOutput', false));
            hold off;
        end
    end
end
