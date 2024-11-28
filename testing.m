LEFT_BOUNDARY = -20;
RIGHT_BOUNDARY = 20;
NUM_OF_POINTS = 1000;
TOT_TIME = 10 * pi;
NUM_OF_STEPS = 1000;

% potential parameter
A = 0.01;
OMEGA = 0.5;

% initialiatState (SHO)
initialStateParams.coefficients = [1]; % complex coefficient of each basis
initialStateParams.bases = [0];       % quantum number for |0> and |1>
initialStateParams.center = 0;           % center of oscillator
initialStateParams.M = 1;                % mass of oscillator
initialStateParams.OMEGA = 1;            % angular freq of oscillator
initialStateParams.HBAR = 1;             % hbar set to 1

targetBasis.coefficients = [1];   % complex coefficient of each basis
targetBasis.bases = [1];          % quantum number for |0> and |1>
targetBasis.center = 0;           % center of oscillator
targetBasis.M = 1;                % mass of oscillator
targetBasis.OMEGA = 1;            % angular freq of oscillator
targetBasis.HBAR = 1;             % hbar set to 1

sim = QuantumDynamicsSimulator(LEFT_BOUNDARY, RIGHT_BOUNDARY, NUM_OF_POINTS, TOT_TIME, NUM_OF_STEPS);
sim = sim.initializeSHOState(initialStateParams);
sim = sim.runSimulationSplitOperatorMethod(A,OMEGA);
sim = sim.plotResults();
sim = sim.initializeSHOState(targetBasis);
[sim, fig1] = sim.plotTransitionProbability(targetBasis);
sim.plotTransition(A, OMEGA, 1, 0, TOT_TIME);
sim = sim.animateEvolution();
