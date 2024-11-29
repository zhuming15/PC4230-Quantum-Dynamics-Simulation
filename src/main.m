%% Main Script: Quantum Dynamics Simulation
% Clear workspace and close all figures
clear; clc; close all;

%% Define Parameters
% Spatial grid
leftBoundary = -20;
rightBoundary = 20;
numPoints = 1000;

% Time parameters
totalTime = 10 * pi;
numSteps = 1000;
frameStep = 1;

% Driving potential parameters
A = 0.01;          % Amplitude of driving potential
OMEGA = 0.5;       % Driving frequency

% Harmonic oscillator parameters
M = 1;             % Mass
OMEGA_HO = 1;      % Harmonic oscillator frequency
HBAR = 1;          % Reduced Planck's constant

%% Initialize GridManager
gridManager = GridManager(leftBoundary, rightBoundary, numPoints);

%% Initialize Initial Wavefunction (Ground State of Harmonic Oscillator)
wavefunction = Wavefunction(gridManager);
wavefunction = wavefunction.initializeHarmonicOscillator(struct( ...
    'bases', 0, ...          % Ground state
    'coefficients', 1, ...
    'center', 0, ...
    'M', M, ...
    'OMEGA', OMEGA_HO, ...
    'HBAR', HBAR ...
    ));

%% Define Propagators
kineticEnergyFunc = @(p) p.^2 / 2; % Kinetic energy
potentialEnergyFunc = @(x, t) 0.5 * x.^2 + A * sin(x) .* cos(OMEGA * t); % Harmonic potential + driving
propagators = Propagators(gridManager, kineticEnergyFunc, potentialEnergyFunc);

%% Create and Run Quantum Simulator
simulator = QuantumSimulator(gridManager, wavefunction, propagators, totalTime, numSteps);
simulator = simulator.runSimulation(frameStep);

%% Define Target Wavefunction (First Excited State of Harmonic Oscillator)
targetWavefunction = Wavefunction(gridManager).initializeHarmonicOscillator(struct( ...
    'bases', 1, ...          % First excited state
    'coefficients', 1, ...
    'center', 0, ...
    'M', M, ...
    'OMEGA', OMEGA_HO, ...
    'HBAR', HBAR ...
    ));

% simulator.animateWavefunction(frameStep);
%% Plot Numerical Transition Probabilities
numericalFig = simulator.plotNumericalTransitionProbability(targetWavefunction);

%% Plot Analytical Transition Probabilities
% Driving potential and spatial grid
V = (A )%* sin(gridManager.getSpatialGrid()));

% Analytical transition probability plot
analyticalFig = simulator.plotAnalyticalTransitionProbability( ...
    V, OMEGA, targetWavefunction);

%% Optional: Animate Wavefunction Evolution
% Uncomment to create an animation
% frameRate = 30; % Frames per second
% savePath = 'wavefunction_evolution.mp4'; % Video file path
% simulator.animateWavefunction(frameRate, savePath);