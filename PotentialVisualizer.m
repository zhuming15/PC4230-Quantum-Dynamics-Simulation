classdef PotentialVisualizer
    properties
        potentialFunc  % Handle to the potential function, a function of (x, t)
        spatialDomain  % Array defining the range of x values
        timeDomain     % Array defining the range of t values for animation
    end
    
    methods
        function obj = PotentialVisualizer(func, xRange, tRange)
            % Constructor to initialize the PotentialVisualizer
            obj.potentialFunc = func;
            obj.spatialDomain = linspace(xRange(1), xRange(2), 400);
            obj.timeDomain = linspace(tRange(1), tRange(2), 100);
        end
        
        function [obj, figHandle] = plotPotential(obj, time)
            % Method to plot the potential function at a specific time
            V = obj.potentialFunc(obj.spatialDomain, time);
            figHandle = figure;
            plot(obj.spatialDomain, V, 'LineWidth', 2);
            title(sprintf('Potential at t = %.2f', time));
            xlabel('Spatial coordinate, x');
            ylabel('Potential V(x, t)');
            grid on;
        end
        
        function [obj, figHandle] = animatePotential(obj)
            % Method to animate the potential function over time
            figHanlde = figure;
            hPlot = plot(obj.spatialDomain, obj.potentialFunc(obj.spatialDomain, obj.timeDomain(1)), 'LineWidth', 2);
            title('Dynamic Potential V(x, t)');
            xlabel('Spatial coordinate, x');
            ylabel('Potential V(x, t)');
            grid on;

            % Ensure ylim is set to accommodate all potential fluctuations over time
            allValues = arrayfun(@(t) obj.potentialFunc(obj.spatialDomain, t), obj.timeDomain, 'UniformOutput', false);
            minY = min(cellfun(@min, allValues)-5);
            maxY = max(cellfun(@max, allValues));
            ylim([minY - 1, maxY + 1]);

            % Update plot for each time step
            for t = obj.timeDomain
                V = obj.potentialFunc(obj.spatialDomain, t);
                set(hPlot, 'YData', V);
                title(sprintf('Dynamic Potential V(x, t) at t = %.2f', t));
                drawnow;
                pause(0.05); % Adjust for desired animation speed
            end
        end
    end
end
