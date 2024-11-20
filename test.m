%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Time Evolution of a Quantum Particle in a 1D Harmonic Trap with
%   Time-Periodic Perturbation using Split-Operator Fourier Method
%   All quantities are in dimensionless units
%   Unit of energy: hbar*omega
%   Unit of length: l = sqrt(hbar/(m*omega))
%   Unit of momentum: p = hbar/l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear Workspace and Initialize Timer
clear all;            % Remove all variables from the workspace
clc;                  % Clear the command window
close all;            % Close all figure windows
tic;                  % Start timer

%% Define Simulation Parameters
a = -15;              % Left end point of the trap
b = +15;              % Right end point of the trap
L = b - a;            % Width of the trap (60 units)
N = 1024;             % Number of spatial points (power of 2 for FFT efficiency)
X = a + L * (0:N-1) / N;  % Dimensionless position coordinates [0 to N-1]
dx = L / N;               % Spatial step size
P = (2 * pi / L) * [0:N/2-1, -N/2:-1];  % Dimensionless momentum grid
dt = 0.1;                % Time step size (adjust for convergence)
T_total = 30;              % Total time duration of the evolution
M = T_total / dt;          % Total number of time steps
% T=5*pi;                         % Time duration of the evolution
% M = 10^3;                     % Total No. of steps in the evolution
% dt = T/M;
%% Define Potential Parameters
A = 0.0;                  % Driving amplitude
omega = 0.5;              % Driving frequency

%% Define Propagators
% Kinetic Propagator (Momentum Space)
UT = exp(-1i * (P.^2 / 2) * dt);   % Kinetic energy propagator

%% Initialize Wavepacket in the Ground State
% For the harmonic oscillator, the ground state is a Gaussian:
sigma = 1.0;              % Width of the Gaussian (matches harmonic oscillator ground state)
X0 = 1;                 % Centered at origin
psi_initial = exp(-(X(1:N)-X0).^2/(2*sigma^2));  % Gaussian wavepacket
psi_initial_norm = psi_initial / sqrt(sum(abs(psi_initial).^2));  % Normalize
psi = psi_initial;     % Initialize wavepacket
%% Prepare Figure for Animation
figure;
h = plot(X, abs(psi).^2, 'b', 'LineWidth', 1.5);  % Plot the initial probability density in blue
xlabel('Position X');
ylabel('Probability Density |\psi(x)|^2');
title('Wavepacket Evolution in a 1D Harmonic Trap with Time-Periodic Perturbation');
ylim([0, max(abs(psi).^2) * 1.2]);               % Set y-axis limits for better visibility
grid on;
hold on;

%% Prepare Video Writer (Optional)
save_animation = true;                            % Set to true to save the animation
if save_animation
    video_filename = 'wavepacket_evolution.mp4';  % Name of the output video file
    v = VideoWriter(video_filename, 'MPEG-4');    % Create a VideoWriter object
    v.FrameRate = 30;                             % Set frame rate
    open(v);                                       % Open the video file for writing
end

%% Time Evolution Loop with Real-Time Plot Updates
for m = 1:M
    current_time = m * dt;  % Current time
    
    % Define Time-Dependent Perturbation V(t) = A sin(x) cos(omega t)
    V_pert = A * sin(X) * cos(omega * current_time);  % Perturbation at current time
    V_total = 0.5 * X.^2 + V_pert;                   % Total potential: harmonic + perturbation
    
    % Potential Propagator for Full Step (Exponentiated Potential)
    UV = exp(-1i * V_total * dt);                   % Potential evolution operator
    
    % Split-Step Method: Apply half potential, full kinetic, half potential
    % Step 1: Half-step potential evolution
    psi = (UV).^0.5 .* psi;                           % Apply sqrt(UV) to represent half-step
    
    % Step 2: Full-step kinetic evolution
    phi = fft(psi);                                  % Fourier transform to momentum space
    phi = UT .* phi;                                  % Apply kinetic propagator
    psi = ifft(phi);                                  % Inverse Fourier transform back to position space
    
    % Step 3: Half-step potential evolution
    psi = (UV).^0.5 .* psi;                           % Apply sqrt(UV) again for the second half-step
    
    % Normalize the wavefunction to prevent numerical drift
    psi = psi / sqrt(sum(abs(psi).^2) * dx);        % Re-normalize
    
    % Update Plot at Every Time Step
    set(h, 'YData', abs(psi).^2);            % Update probability density
    drawnow;                                   % Render the updated plot immediately
    
    % Capture Frame for Video (Optional)
    if save_animation
        frame = getframe(gcf);                 % Capture the current figure as a frame
        writeVideo(v, frame);                  % Write the frame to the video file
    end
end

%% Finalize Video (If Recording)
if save_animation
    close(v);                                       % Close the video writer object
    disp(['Animation saved as ', video_filename]);  % Display confirmation
end

%% Plot Final State
figure;
plot(X, abs(psi).^2, 'r', 'LineWidth', 1.5);      % Plot the final probability density in red
xlabel('Position X');
ylabel('Probability Density |\psi|^2');
title('Final Wavepacket After Evolution');
ylim([0, max(abs(psi).^2) * 1.2]);               % Set y-axis limits
grid on;