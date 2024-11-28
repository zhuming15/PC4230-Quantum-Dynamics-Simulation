classdef Propagators
    % PROPAGATORS Handles kinetic and potential propagators for quantum systems.

    properties
        kineticEnergyFunc     % Function handle for kinetic energy
        potentialEnergyFunc   % Function handle for potential energy
    end

    methods
        function obj = Propagators(kineticEnergyFunc, potentialEnergyFunc)
            obj.kineticEnergyFunc = kineticEnergyFunc;
            obj.potentialEnergyFunc = potentialEnergyFunc;
        end

        function UT = generateKineticPropagator(obj, momentumGrid, dt)
            % Generate the kinetic propagator
            UT = exp(-1i * obj.kineticEnergyFunc(momentumGrid) * dt);
        end

        function UV_list = generatePotentialPropagators(obj, spatialGrid, dt, totalTime, numSteps)
            % Generate the potential propagators for each time step
            times = linspace(0, totalTime, numSteps);
            UV_list = cell(1, numSteps); % Preallocate cell array
            for k = 1:numSteps
                potentialEnergy = obj.potentialEnergyFunc(spatialGrid, times(k));
                UV_list{k} = exp(-1i * potentialEnergy * dt / 2);
            end
        end
    end
end
