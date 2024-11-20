% main.m - Script to use the WavePacketSimulation class

clear;
clc;

% Step 1: Create an instance of WavePacketSimulation
% Parameters: a, b, N, T, M
sim = QuantumDynamicsSimulator(-80, 80, 10000, 4*pi, 1000);
sim = sim.initializeWavePacket(20, 1)

% Step 2: Run the simulation
sim = sim.runSimulation(6,0.5);

% Step 3: Plot the initial and final states
% [sim, fig] = sim.plotResults();

% Step 4: Animate the wavepacket evolution
[sim, animFig] = sim.animateWavePacket();
