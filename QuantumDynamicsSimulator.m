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
            obj.stateEvolution = complex(zeros(points, steps));
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
        
        function [basis, energy, obj] = generateHermiteBasis(obj, n, center, M, OMEGA, HBAR)        
            x = obj.spatialGrid(1:obj.numPoints) - center;
            H_n = hermiteH(n, sqrt(M * OMEGA / HBAR) * x);
            basis = H_n .* exp(-M * OMEGA * x.^2 / (2 * HBAR));
            basis = basis / sqrt(sum(abs(basis).^2));
            energy = HBAR * OMEGA * (n + 0.5);
        end

        function obj = initializeSHOState(obj, params)
            obj.initState = zeros(size(obj.spatialGrid)); 
            for k = 1:length(params.bases)
                n = params.bases(k);
                coeff = params.coefficients(k);
                [basis, ~] = obj.generateHermiteBasis(n, params.center, params.M, params.OMEGA, params.HBAR);
                obj.initState = obj.initState + coeff * basis;
            end
            obj.initState = obj.initState / sqrt(sum(abs(obj.initState).^2));
        end

        function obj = runSimulationSplitOperatorMethod(obj, A, omega)
            % Pre-compute exponential terms for the kinetic energy part
            % Kinetic term Propagator in momentum space
            UT = exp(-1i * (obj.momentumGrid.^2 / 2) * obj.dt);

            % Potential term in real space
            times = linspace(0, obj.totalTime, obj.numSteps);
            UV_list = arrayfun(@(t) exp(-1i * ((obj.spatialGrid.^2 / 2) + ...
                A * sin(obj.spatialGrid) * cos(omega * t)) * obj.dt / 2), times, 'UniformOutput', false);

            cur_psi = obj.initState;
            obj.stateEvolution(:, 1) = cur_psi;  % Store the state at this time step
            for m = 2:obj.numSteps
                UV = UV_list{m};
                cur_psi = UV .* cur_psi;  % Apply half-step potential in real space
                cur_psi = fft(cur_psi);  % Transform to momentum space
                cur_psi = UT .* cur_psi;  % Apply full-step kinetic in momentum space
                cur_psi = ifft(cur_psi);  % Transform back to real space
                cur_psi = UV .* cur_psi;  % Apply second half-step potential in real space
        
                obj.stateEvolution(:, m) = cur_psi;  % Store the state at this time step
            end
            obj.finalState = cur_psi;
        end


        function [obj, figHandle] = plotTransitionProbability(obj, targetBasis)
            % basis
            basis = obj.generateHermiteBasis( ...
                targetBasis.bases, ...
                targetBasis.center, ...
                targetBasis.M, ...
                targetBasis.OMEGA, ...
                targetBasis.HBAR ...
                );
        
            % Transition probability list
            transitionProbabilities = zeros(1,obj.numSteps);
            dx = obj.spatialGrid(2) - obj.spatialGrid(1);

            % Compute probabilities for each timestep
            for t = 1:obj.numSteps
                stateAtT = obj.stateEvolution(:, t);
                transitionProbabilities(t) = abs(sum(conj(basis) * stateAtT))^2;
            end

            % Create time vector
            time = linspace(0, obj.totalTime, obj.numSteps);

            % Plot transition probabilities
            figHandle = figure;
            plot(time, transitionProbabilities, 'LineWidth', 2);

            % Customize plot
            title('Transition Probability vs. Time');
            xlabel('Time');
            ylabel('Transition Probability');
            % ylim([0,1]);
            legend show;
            grid on;
        end

        function obj = plotResults(obj)
            % Plot initial and final probability densities
            figure;
            subplot(2, 1, 1);
            plot(obj.spatialGrid, obj.computeProbDist(obj.initState), 'LineWidth', 2);
            title('Initial State Probability Density');
            xlabel('Position'); ylabel('Density'); grid on;

            subplot(2, 1, 2);
            plot(obj.spatialGrid, obj.computeProbDist(obj.finalState), 'LineWidth', 2);
            title('Final State Probability Density');
            xlabel('Position'); ylabel('Density'); grid on;
        end
        
        function obj = animateEvolution(obj)
            % Animate wavepacket evolution
            figure;
            hPlot = plot(obj.spatialGrid, obj.computeProbDist(obj.stateEvolution(:, 1)), 'LineWidth', 2);
            ylim([0, max(abs(obj.stateEvolution(:)).^2) * 1.2]);
            xlabel('Position'); ylabel('Probability Density'); title('Wavepacket Evolution');
            for m = 2:obj.numSteps
                set(hPlot, 'YData', obj.computeProbDist(obj.stateEvolution(:, m)));
                drawnow;
            end
        end

        function [obj, figHandle] = plot3DStateEvolution(obj)
            figHandle = figure;
            [X, T] = meshgrid(obj.spatialGrid, linspace(0, obj.totalTime, obj.numSteps));
            surf(X, T, abs(obj.stateEvolution').^2, 'EdgeColor', 'none');
            title('Wavepacket Evolution');
            xlabel('Position');
            ylabel('Time');
            zlabel('Probability Amplitude');
            view(3);
        end

        function [obj, figHandle] = plotTransition(obj, A, omega, omega_nm, timeStart, timeEnd)
            % Plot transition probability vs time for potential sin(x)
            potFunc = A * sin(obj.spatialGrid);
            timeVector = linspace(timeStart, timeEnd, obj.numPoints);
            transProbCurve = zeros(1, obj.numPoints);
            
            for i = 1:obj.numPoints
                t = timeVector(i);
                transProbCurve(i) = obj.computeAnalyTransProb(potFunc, omega, omega_nm, t);
            end

            figHandle = figure;
            plot(timeVector, transProbCurve, 'LineWidth', 2);
            title('Transition Probability vs Time');
            xlabel('Time'); ylabel('Analytic Transition Probability'); grid on;
        end

        function transProb = computeAnalyTransProb(obj, potFunc, omega, omega_nm, t)
            % Compute analytical transition probability
            potentialTerm = 4 * abs(sum(potFunc))^2;
            freqDiff = omega - omega_nm;
            freqFactor = (sin(freqDiff * t / 2) / freqDiff)^2;
            transProb = potentialTerm * freqFactor;
        end

        function [probDist, obj] = computeProbDist(obj, ket1)
            probDist = abs(ket1).^2;
        end

        function [transProb, obj] = computeTransProb(obj, ket1, ket2)
            transProb = abs(sum(conj(ket1) .* ket2))^2;
        end
    end
end