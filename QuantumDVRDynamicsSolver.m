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
            obj.spatialGrid = linspace(left, right, numPts-1);
            obj.hamiltonian = [];  % Initialize to empty; will be calculated later
            obj.eigenvalues = [];
            obj.eigenfunctions = [];
        end
        
        function obj = generateHamiltonian(obj, potentialFunc)
            % Generates the Hamiltonian matrix using the DVR approach and specified potential
            disp(obj.spatialGrid);
            obj.hamiltonian = zeros(obj.numPoints-1);
            kinetic = zeros(obj.numPoints-1);
            potential = diag(arrayfun(potentialFunc, obj.spatialGrid));

            c = (1/2)*(pi/(obj.rightBoundary - obj.leftBoundary))^2*(2/obj.numPoints);
            nsum = 1:obj.numPoints-1;
            nsquare = (1:obj.numPoints-1).^2;
            % Assembling the Hamiltonian matrix
            for i = 1:obj.numPoints-1
                for j = i:obj.numPoints-1
                    kinetic(i, j) = c*sum(nsquare.*sin(nsum*pi*i/obj.numPoints).*sin(nsum*pi*j/obj.numPoints));
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
            title('Eigenfunctions of the quantum system');
            xlabel('Spatial coordinate');
            ylabel('Wavefunction');
            legend(arrayfun(@(x) ['\psi_', num2str(x)], EigenfunctionsToPlot, 'UniformOutput', false));
            hold off;
        end

        function [obj, figHandle] = plotEigenenergies(obj, eigenIndices)
            % Plots eigenenergies (eigenvalues) for specified indices
            if isempty(obj.eigenvalues)
                error('Eigenvalues have not been calculated. Call solveEigenproblems first.');
            end
            if max(eigenIndices) > length(obj.eigenvalues)
                error('Index exceeds the number of calculated eigenvalues.');
            end

            figHandle = figure;  % Create a new figure window
            energies = obj.eigenvalues(eigenIndices);  % Extract the eigenvalues for given indices
            disp(energies);
            plot(eigenIndices, energies, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
            title('Eigenenergies of the Quantum System');
            xlabel('Eigenvalue Index');
            ylabel('Quantum Energy Levels');  % Updated label for clarity
            grid on;
            % Set the axes to have ticks at each index for clarity
            ax = gca;
            ax.XTick = eigenIndices;
            ax.XLim = [min(eigenIndices) - 1, max(eigenIndices) + 1];
        end
    end
end
