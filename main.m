% main.m - Script to use the WavePacketSimulation class

clear;
clc;

% runSimulation()
% solveEigenvalEigenfunc()
% plotPotentialGraph()


% Methodology
% 1. Find Eigenenergies of Harmonic Oscillator using Hermite Polynomials (Theoretical Value)
% 2. Find Eigenenergies of Harmonic Oscillator using DVR (Numerical Value)
% 3. Compure and ensure DVR representation is accurate
% 4. Prepare the initial state using DVR representation, ensure the parameter is standarized.
% 5. Implement function to compute transition probability. 
% 6. Implement function to plot/animate the transition probability.
% 7. Consider adiabatic limit of omega -> 0, and compare the transition probability with the theoretical value.
% 8. Consider omega = 0, A -> infinity and compare with mulitple finite square well.
% 9. Consider omega -> inifnity. State should not be able to reponse to the time dependent potential.

% TODO - Tidy up the plot, include A and omega label in plot

% Define the potential functions
function plotPotentialGraph
    % Plot parameter
    LEFT_BOUNDARY = -50
    RIGHT_BOUNDARY = 50
    START_TIME = 0
    END_TIME = 10 * pi

    % Potential parameter
    A = 100
    OMEGA = 0.5
    
    harmonicPotential = @(x) x.^2 / 2;  % Correct, assumes x is an array
    timeDependentPotential = @(x, A, omega, t) A .* sin(x) .* cos(omega * t);
    potentialFunc = @(x, t) harmonicPotential(x) + timeDependentPotential(x, A, OMEGA, t);

    visualizer = PotentialVisualizer( ...
        potentialFunc, ...
        [LEFT_BOUNDARY, RIGHT_BOUNDARY], ...
        [START_TIME, END_TIME] ...
        );
    visualizer.animatePotential();
end

function solveEigenvalEigenfunc
    % DVR parameter
    LEFT_BOUNDARY = -20
    RIGHT_BOUNDARY = 20
    NUM_OF_POINTS = 500

    % Potential parameter
    A = 0
    OMEGA = 2
    TIME = 2 *pi

    % Plot parameter
    EIGENFUCNTION_TO_PLOTS = [1, 2, 3, 4, 5];

    % Initialize the QuantumDVRDynamicsSolver with DVR parameter
    DVRSolver = QuantumDVRDynamicsSolver(LEFT_BOUNDARY, RIGHT_BOUNDARY, NUM_OF_POINTS);
    
    % Generate the Hamiltonian matrix
    harmonicPotential = @(x) x.^2 / 2;
    timeDependentPotential = @(x, A, omega, t) A .* sin(x) .* cos(omega * t);
    potentialFunc = @(x) harmonicPotential(x) + timeDependentPotential(x, A, OMEGA, TIME);
    DVRSolver = DVRSolver.generateHamiltonian(potentialFunc);
    
    % Solve the eigenvalue problem to find eigenstates and eigenvalues
    [DVRSolver, vec, eneg] = DVRSolver.solveEigenproblems();
    
    % Plot the eigenfunctions and eigenvalues
    [DVRSolver, fig] = DVRSolver.plotEigenfunctions(EIGENFUCNTION_TO_PLOTS);
    [DVRSolver, fig] = DVRSolver.plotEigenenergies(EIGENFUCNTION_TO_PLOTS);
end


function runSimulation
    % simulation parameter
    LEFT_BOUNDARY = -50
    RIGHT_BOUNDARY = 50
    NUM_OF_POINTS = 1000
    TOT_TIME = 10 * pi
    NUM_OF_STEPS = 10000

    %initial state parameter of Guassian Wavapacket
    X_0 = 0
    SIGMA = 1

    % potential parameter
    A = 1000
    OMEGA = 0.5
    
    sim = QuantumDynamicsSimulator(LEFT_BOUNDARY, RIGHT_BOUNDARY, NUM_OF_POINTS, TOT_TIME, NUM_OF_STEPS);
    sim = sim.initializeWavePacket(X_0, SIGMA)
    sim = sim.runSimulationSplitOperatorMethod(A,OMEGA);
    
    [sim, fig] = sim.plotResults();
    
    [sim, animFig] = sim.animateWavePacket();
end
 