classdef TransitionAnalysis
    % TRANSITIONANALYSIS Handles transition probability calculations for
    % numerical and analytical methods.

    methods (Static)
        function prob = computeNumericalTransition(stateAtT, targetState)
            prob = abs((conj(targetState(:))' * stateAtT)).^2;
        end

        % function prob = computeAnalyticalTransition(potential, omega, omega_nm, t)
        % 
        %     % Compute potential term (constant factor)
        %     potentialTerm = 4 * abs(sum(potential))^2;
        % 
        %     % Frequency difference factor
        %     freqDiff = omega - omega_nm;
        %     if abs(freqDiff) < 1e-10
        %         freqFactor = (t / 2)^2; % Special case for resonance (freqDiff ~ 0)
        %     else
        %         freqFactor = (sin(freqDiff * t / 2) / freqDiff)^2;
        %     end
        % 
        %     % Compute transition probability
        %     prob = potentialTerm * freqFactor;
        % end

        function fig = plotNumericalTransitionProbability(stateEvolution, targetState, timeVector)
            transitionProbs = TransitionAnalysis.computeNumericalTransition(stateEvolution, targetState);

            fig = figure;
            plot(timeVector, transitionProbs, 'LineWidth', 2);
            title('Transition Probability (Numerical)');
            xlabel('Time'); ylabel('Transition Probability'); grid on;
        end

        % function fig = plotAnalyticalTransition(potential, omega, omega_nm, spatialGrid, timeVector)
        %     % Compute transition probabilities over all time points
        %     transitionProbs = arrayfun(@(t) TransitionAnalysis.computeAnalyticalTransition( ...
        %         potential, omega, omega_nm, t), timeVector);
        % 
        %     % Plot transition probabilities
        %     fig = figure;
        %     plot(timeVector, transitionProbs, 'LineWidth', 2);
        %     title('Transition Probability (Analytical)');
        %     xlabel('Time'); ylabel('Transition Probability'); grid on;
        % end
    end
end
