classdef QuantumSimulator
    % QUANTUMSIMULATOR Main simulation engine for quantum dynamics.

    properties
        gridManager      % Instance of GridManager
        wavefunction     % Instance of Wavefunction
        propagator       % Instance of Propagators
        stateEvolution   % Matrix storing wavefunctions at selected frames (state is a column vector)
        totalTime        % Total simulation time
        dt               % Time step
        numSteps         % Number of time steps
    end

    methods
        %% Constructor
        function obj = QuantumSimulator(gridManager, wavefunction, propagator, totalTime, numSteps)
            % Constructor to initialize the simulation engine
            obj.gridManager = gridManager;
            obj.wavefunction = wavefunction;
            obj.propagator = propagator;
            obj.totalTime = totalTime;
            obj.numSteps = numSteps;
            obj.dt = totalTime / numSteps;
        end

        function obj = runSimulation(obj, frameStep)
            % Run the split-operator simulation with selective frame storage
            kineticProp = obj.propagator.generateKineticPropagator(obj.dt);
            potentialPropList = obj.propagator.generatePotentialPropagators( ...
                obj.dt, obj.totalTime, obj.numSteps);
        
            % Preallocate memory for selected frames
            numFrames = ceil(obj.numSteps / frameStep);
            obj.stateEvolution = complex(zeros(obj.gridManager.getNumPoints(), numFrames));
        
            % Initialize state evolution
            currentState = obj.wavefunction.state;
            frameIndex = 1;
            obj.stateEvolution(:, frameIndex) = currentState;

            for step = 2:obj.numSteps
                currentState = obj.splitOperatorStep(currentState, kineticProp, potentialPropList{step});

                % Store wavefunction at specified frame intervals
                if mod(step, frameStep) == 0
                    frameIndex = frameIndex + 1;
                    obj.stateEvolution(:, frameIndex) = currentState;
                end
            end
            obj.wavefunction.state = currentState;
        end

        function fig = plotNumericalTransitionProbability(obj, targetState)
            numFrames = size(obj.stateEvolution, 2);
            timeVector = linspace(0, obj.totalTime, obj.numSteps);
            fig = TransitionAnalysis.plotNumericalTransitionProbability( ...
                obj.stateEvolution, ...
                targetState.state, ...
                timeVector ...
                );
        end

        function fig = plotAnalyticalTransitionProbability(obj, V_nm, omega, targetState)
            omega_nm = abs(targetState.params.bases - obj.wavefunction.params.bases);
            disp(omega_nm);
            numFrames = size(obj.stateEvolution, 2)
            timeVector = linspace(0, obj.totalTime, numFrames);
            fig = TransitionAnalysis.plotAnalyticalTransitionProbability( ...
                V_nm, ...
                omega, ...
                omega_nm, ...
                timeVector ...
                );
        end

        function animateWavefunction(obj, frameRate, savePath)
            % ANIMATEWAVEFUNCTION Animates the time evolution of the wavefunction
            % with the potential curve overlaid.
            % INPUTS:
            %   frameRate - Frame rate for the animation (frames per second)
            %   savePath - (Optional) Path to save the animation as a video file
            %              If omitted, the animation will not be saved.

            % Extract spatial grid and frames
            spatialGrid = obj.gridManager.getSpatialGrid();
            numFrames = size(obj.stateEvolution, 2);

            % Generate potential curve for overlay
            potentialFunc = obj.propagator.potentialEnergyFunc; % Potential function
            timeGrid = linspace(0, obj.totalTime, numFrames);
            potentialOverlay = potentialFunc(spatialGrid, 0); % Time-independent part for static overlay

            % Set up the figure
            fig = figure('Name', 'Wavefunction Evolution', 'NumberTitle', 'off');
            fig.Position = [100, 100, 1200, 400];

            % Subplots: real/imaginary part with potential overlay and probability density
            ax1 = subplot(1, 2, 1);
            hold(ax1, 'on');
            realPartPlot = plot(ax1, spatialGrid, real(obj.stateEvolution(:, 1)), 'b', 'LineWidth', 1.5, 'DisplayName', 'Real Part');
            imagPartPlot = plot(ax1, spatialGrid, imag(obj.stateEvolution(:, 1)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Imaginary Part');
            potentialPlot1 = plot(ax1, spatialGrid, potentialOverlay, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Potential');
            title(ax1, 'Real/Imaginary Parts with Potential Overlay');
            xlabel(ax1, 'Position');
            ylabel(ax1, 'Amplitude / Potential');
            legend(ax1, 'Location', 'best');
            grid(ax1, 'on');
            ylim(ax1, [-1, 1] * max(abs(obj.stateEvolution(:)))); % Dynamic range

            ax2 = subplot(1, 2, 2);
            hold(ax2, 'on');
            probDensityPlot = plot(ax2, spatialGrid, abs(obj.stateEvolution(:, 1))^2, 'r', 'LineWidth', 1.5);
            potentialPlot2 = plot(ax2, spatialGrid, potentialOverlay, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Potential');
            title(ax2, 'Probability Density');
            xlabel(ax2, 'Position');
            ylabel(ax2, 'Density');
            grid(ax2, 'on');
            ylim(ax2, [0, max(abs(obj.stateEvolution(:)).^2)]); % Dynamic range

            % Video writer (if savePath is specified)
            if nargin > 2
                videoWriter = VideoWriter(savePath, 'MPEG-4');
                videoWriter.FrameRate = frameRate;
                open(videoWriter);
            end

            % Animation loop
            for frame = 1:numFrames
                % Update potential overlay for the current time (if time-dependent)
                currentTime = timeGrid(frame);
                potentialOverlay = potentialFunc(spatialGrid, currentTime); % Time-dependent potential
                potentialPlot1.YData = potentialOverlay;
                potentialPlot2.YData = potentialOverlay;

                % Update wavefunction plots
                realPartPlot.YData = real(obj.stateEvolution(:, frame));
                imagPartPlot.YData = imag(obj.stateEvolution(:, frame));
                probDensityPlot.YData = abs(obj.stateEvolution(:, frame)).^2;


                % Update subplot titles with time information
                ax1.Title.String = sprintf('Real/Imaginary Parts (t = %.2f)', currentTime);
                ax2.Title.String = sprintf('Probability Density (t = %.2f)', currentTime);

                % Draw and capture the frame
                drawnow;

                % Write the frame to video if needed
                if nargin > 2
                    frameData = getframe(fig);
                    writeVideo(videoWriter, frameData);
                end
            end

            % Close video writer if used
            if exist('videoWriter', 'var')
                close(videoWriter);
            end
        end
    end

    methods (Access = private)
        function updatedState = splitOperatorStep(~, state, kineticProp, potentialProp)
            updatedState = potentialProp .* state;           % Apply half potential
            updatedState = fft(updatedState);                % Transform to momentum space
            updatedState = kineticProp .* updatedState;      % Apply kinetic propagator
            updatedState = ifft(updatedState);               % Transform back to real space
            updatedState = potentialProp .* updatedState;    % Apply second half potential
        end
    end
end
