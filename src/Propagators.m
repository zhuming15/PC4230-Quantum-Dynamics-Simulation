classdef Propagators
    % PROPAGATORS Handles kinetic and potential propagators for quantum systems.

    properties
        gridManager          % Instance of GridManager
        kineticEnergyFunc    % Function handle for kinetic energy
        potentialEnergyFunc  % Function handle for potential energy
    end

    methods
        function obj = Propagators(gridManager, kineticEnergyFunc, potentialEnergyFunc)
            % Constructor to initialize the propagators
            obj.gridManager = gridManager;
            obj.kineticEnergyFunc = kineticEnergyFunc;
            obj.potentialEnergyFunc = potentialEnergyFunc;
        end

        function kineticProp = generateKineticPropagator(obj, dt)
            % Generate the kinetic propagator
            pGrid = obj.gridManager.getMomentumGrid(); % Access via getter
            kineticProp = exp(-1i * obj.kineticEnergyFunc(pGrid) * dt);
        end

        function potentialPropList = generatePotentialPropagators(obj, dt, totalTime, numSteps)
            % Generate the potential propagators for the entire simulation
            spatialGrid = obj.gridManager.getSpatialGrid(); % Access via getter
            timeGrid = linspace(0, totalTime, numSteps);

            % Preallocate cell array for propagators
            potentialPropList = cell(1, numSteps);

            for k = 1:numSteps
                potentialEnergy = obj.potentialEnergyFunc(spatialGrid, timeGrid(k));
                potentialPropList{k} = exp(-1i * potentialEnergy * dt / 2);
            end
        end
    end
end
