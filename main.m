% main.m - Script to use the WavePacketSimulation class

clear;
clc;

% runSimulation()
solveEigenvalEigenfunc()
% plotPotentialGraph()

% TODO - Tidy up the plot, include A and omega label in plot
% Define the potential functions
function plotPotentialGraph
    % Plot parameter
    LEFT_BOUNDARY = -50
    RIGHT_BOUNDARY = 50
    START_TIME = 0
    END_TIME = 2 * pi

    % Potential parameter
    A = 6
    OMEGA = 1
    
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
    OMEGA = 1
    TIME = pi / 2

    % Plot parameter
    EIGENFUCNTION_TO_PLOTS = [1,2,3];

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
    TOT_TIME = 2 * pi
    NUM_OF_STEPS = 1000

    %initial state parameter of Guassian Wavapacket
    X_0 = 10
    SIGMA = 2

    % potential parameter
    A = 0
    OMEGA = 1
    
    sim = QuantumDynamicsSimulator(LEFT_BOUNDARY, RIGHT_BOUNDARY, NUM_OF_POINTS, TOT_TIME, NUM_OF_STEPS);
    sim = sim.initializeWavePacket(X_0, SIGMA)
    sim = sim.runSimulationSplitOperatorMethod(A,OMEGA);
    
    [sim, fig] = sim.plotResults();
    
    [sim, animFig] = sim.animateWavePacket();
end
 