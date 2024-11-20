classdef QuantumDynamicsSimulator
    properties
        leftBoundary    % Start of spatial domain
        rightBoundary   % End of spatial domain
        numPoints       % Number of spatial grid points
        totalTime       % Total simulation time
        numSteps        % Number of time steps
        spatialGrid     % Spatial grid array
        momentumGrid    % Momentum grid array
        dt              % Time step size
        initState       % Initial state (wavepacket state)
        finalState      % Final state after simulation
        stateEvolution  % List of state over time
    end
    
    methods
        function obj = QuantumDynamicsSimulator(start, endp, points, time, steps)
            % Constructor to initialize simulation parameters and preallocate space.
            obj.leftBoundary = start;
            obj.rightBoundary = endp;
            obj.numPoints = points;
            obj.totalTime = time;
            obj.numSteps = steps;
            obj.dt = time / steps;
            [obj.spatialGrid, obj.momentumGrid] = obj.generateGrids();
            obj.stateEvolution = complex(zeros(points, steps)); % Initialize complex array
        end
        
        function [X, P] = generateGrids(obj)
            % Generates spatial and momentum grids.
            L = obj.rightBoundary - obj.leftBoundary;
            X = obj.leftBoundary + L * (0:obj.numPoints-1) / obj.numPoints;
            P = (2 * pi / L) * [0:obj.numPoints/2-1, -obj.numPoints/2:-1];
        end
        
        function obj = initializeWavePacket(obj, center, width)
            % Initializes the Gaussian wavepacket.
            prep = exp(-(obj.spatialGrid(1:obj.numPoints) - center).^2 / (2 * width^2));
            obj.initState = prep / sqrt(sum(abs(prep).^2));
        end

        function obj = runSimulationSplitOperatorMethod(obj, A, omega)
            % Pre-compute exponential terms for the kinetic energy part
            % Kinetic term Propagator in momentum space
            UT = exp(-1i * (obj.momentumGrid.^2 / 2) * obj.dt);

            % Potential term in real space
            V = @(A, omega, t) (obj.spatialGrid.^2 / 2) + (A * sin(obj.spatialGrid) * cos(omega * t));
        
            cur_psi = obj.initState;
            for m = 1:obj.numSteps
                % Compute Potential term Propagator in real space
                UV = exp(-1i * V(A, omega, m) * obj.dt / 2);

                cur_psi = UV .* cur_psi;  % Apply half-step potential in real space
                cur_psi = fft(cur_psi);  % Transform to momentum space
                cur_psi = UT .* cur_psi;  % Apply full-step kinetic in momentum space
                cur_psi = ifft(cur_psi);  % Transform back to real space
                cur_psi = UV .* cur_psi;  % Apply second half-step potential in real space
        
                obj.stateEvolution(:, m) = cur_psi;  % Store the state at this time step
            end
            obj.finalState = cur_psi;
        end

        function [obj, figHandle] = plotResults(obj)
            % Plots the initial and final states' probability densities.
            figHandle = figure;
            subplot(2, 1, 1);
            plot(obj.spatialGrid, obj.getProbabilityAmplitude(obj.initState), 'LineWidth', 2);
            title('Initial State Density');
            xlabel('Position');
            ylabel('Density');
            
            subplot(2, 1, 2);
            plot(obj.spatialGrid, obj.getProbabilityAmplitude(obj.finalState), 'LineWidth', 2);
            title('Final State Density');
            xlabel('Position');
            ylabel('Density');
        end
        
        function [obj, figHandle] = animateWavePacket(obj)
            figHandle = figure;
            hPlot = plot(obj.spatialGrid, obj.getProbabilityAmplitude(obj.stateEvolution(:, 1)), 'LineWidth', 2);
            maxY = max(abs(obj.stateEvolution(:)).^2) * 1.2;
            ylim([0, maxY]);
            xlabel('Position');
            ylabel('Probability Amplitute');
            title('Wavepacket Evolution');
            for m = 2:obj.numSteps
                set(hPlot, 'YData', obj.getProbabilityAmplitude(obj.stateEvolution(:, m)));
                drawnow;
            end
            return;
        end
        
        function probAmplitude = getProbabilityAmplitude(obj, wave)
            % Calculate the probability amplitude of a wavepacket
            probAmplitude = abs(wave).^2;
        end
    end
end