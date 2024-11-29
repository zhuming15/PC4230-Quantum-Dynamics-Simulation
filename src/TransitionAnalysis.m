classdef TransitionAnalysis
    % TRANSITIONANALYSIS Handles transition probability calculations for
    % numerical and analytical methods.

    methods (Static)
        function prob = computeNumericalP_nm(stateAtT, targetState)
            prob = abs((conj(targetState(:))' * stateAtT)).^2;
        end
        
        function prob = computeAnalyticalP_nmWithRWA(V_nm, omega, omega_nm, t)
            freqDiff = omega - omega_nm;
            freqFactor = sin(freqDiff * t / 2) / freqDiff;
            prob = 4 * (abs(V_nm) * freqFactor)^2;
        end

        function prob = computeAnalyticalP_nm(V_nm, omega, omega_nm, t)
            term1 = (exp(1i * (omega_nm + omega) * t) - 1) / (omega_nm + omega);
            term2 = (exp(1i * (omega_nm - omega) * t) - 1) / (omega_nm - omega);
            prob = abs(V_nm * (term1 + term2))^2;
        end

        function fig = plotNumericalTransitionProbability(stateEvolution, targetState, timeVector)
            transitionProbs = TransitionAnalysis.computeNumericalP_nm(stateEvolution, targetState);

            fig = figure;
            plot(timeVector, transitionProbs, 'LineWidth', 2);
            title('Transition Probability (Numerical)');
            xlabel('Time'); ylabel('Transition Probability'); grid on;
        end

        function fig = plotAnalyticalTransitionProbability(V_nm, omega, omega_nm, timeVector)
            % Compute transition probabilities over all time points
            transitionProbs = arrayfun(@(t) TransitionAnalysis.computeAnalyticalP_nm( ...
                V_nm, omega, omega_nm, t), timeVector);
            fig = figure;
            plot(timeVector, transitionProbs, 'LineWidth', 2);
            title('Transition Probability (Analytical)');
            xlabel('Time'); ylabel('Transition Probability'); grid on;
        end
    end
end
