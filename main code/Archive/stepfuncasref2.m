%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step Function as Reference System
% Num Met 4 Phys - Homework 7.8.2 (Modified for Step Function)
% MATLAB Code
% Author: Awadhoot Loharkar
% Date: [Date]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization and adding necessary paths
addpath("../../start-kit-student")
mystartdefaults;

% Define constants and parameters
hbar = 1.0545718e-34;    % Reduced Planck's constant (J.s)
me = 9.10938356e-31;     % Electron mass (kg)
qel = 1.60218e-19;       % Electron charge (C)
tol = 1E-12;             % Tolerance for numerical calculations
recipunit = 1E10;        % Reciprocal unit in Å^-1
ekinscale = (hbar * recipunit)^2 / (2 * me * qel); % Conversion factor (Å^-2 to eV)

% Discretization parameters
xmin = -20;   % Start of the simulation domain [Å]
xmax = 100;   % End of the simulation domain [Å]
step = 0.5;   % Spatial step size [Å]
x = xmin:step:xmax;  % Spatial grid
N = length(x);       % Number of discretization points

%% Step 1: Define the Potential Profile (Smoothed Step Function)
% Define the step function U0(x < 0) = 0 and U0(x > 0) = U1 = -0.1 eV
U1 = -0.1;  % Step potential in eV
U0 = zeros(size(x)); % Initialize potential

% Define the smoothed potential profile across the heterostructure
% Modify the potential profile using a hyperbolic tangent function
x_min = 0; 
x_max = 15;
E_bias = abs(U1) / (x_max - x_min); % Electric field associated with bias
U = U0; % Initialize the final potential
U(x >= x_min & x <= x_max) = U1 * (1 - 0.5 * (tanh((x(x >= x_min & x <= x_max) - x_min) / 2) + 1));

% Plot the potential profile
figure;
plot(x, U, 'k', 'LineWidth', 2);
xlabel('Position [Å]');
ylabel('Potential [eV]');
title('Potential Profile with Smoothed Step Function');
grid on;

%% Step 2: Define the Energy Range for Computation
E_min = 0;     % Minimum energy [eV]
E_max = 0.3;   % Maximum energy [eV]
dE = 0.0005;   % Energy step size [eV]
E_range = E_min:dE:E_max;  % Energy range vector

%% Step 3: Compute Reflection, Transmission, and Absorption Coefficients
k = sqrt(2 * me * E_range * qel) / (hbar * recipunit); % Wavevector for each energy
T = zeros(size(E_range)); % Transmission coefficient
R = zeros(size(E_range)); % Reflection coefficient

% Construct Green's function and calculate coefficients
for i = 1:length(E_range)
    k0 = sqrt(E_range(i) / ekinscale);
    G = zeros(N, N);
    % Calculate the Green's function G1 for the smoothed step function potential
    for m = 1:N
        for n = 1:N
            G(m, n) = exp(1i * k0 * abs(x(m) - x(n))) / (2 * 1i * k0);
        end
    end
    
    % Potential weighted matrix
    Waa = diag(U) * step / ekinscale;
    
    % Scattering matrix setup
    ScatMatrix = eye(N) - G * Waa;
    Phi0 = exp(1i * k0 * x); % Incident plane wave
    Phi = ScatMatrix \ transpose(Phi0); % Solve for the scattering wavefunction
    
    % Calculate Reflection and Transmission coefficients
    R(i) = abs((Phi(1) - exp(1i * k0 * x(1))) / exp(1i * k0 * x(1)))^2;
    T(i) = abs(Phi(end) / exp(1i * k0 * x(end)))^2;
end

% Plot Reflection and Transmission coefficients
figure;
plot(E_range, T, 'r', E_range, R, 'b');
xlabel('Energy [eV]');
ylabel('Coefficient');
legend('Transmission', 'Reflection');
title('Transmission and Reflection Coefficients');
grid on;

%% Step 4: Density of States Calculation
lifetime = 1e-9; % Recombination time (lifetime) parameter [s]
gamma = hbar / lifetime; % Broadening factor

% Calculate the density of states using Lorentzian broadening
Density_of_States = (1 / (pi * hbar)) * (gamma ./ (E_range.^2 + gamma^2));

% Plot Density of States
figure;
plot(E_range, Density_of_States, 'g');
xlabel('Energy [eV]');
ylabel('Density of States [1/eV]');
title('Density of States for the Smoothed Step Function Potential');
grid on;

%% Step 5: Tuning the Bias Potential
U1_values = linspace(-0.2, 0.2, 100); % Tuning bias range
T_U1 = zeros(size(U1_values));
R_U1 = zeros(size(U1_values));

for j = 1:length(U1_values)
    U_tuned = U0; % Reinitialize potential
    U_tuned(x >= x_min & x <= x_max) = U1_values(j) * (1 - 0.5 * (tanh((x(x >= x_min & x <= x_max) - x_min) / 2) + 1));
    
    % Recalculate Reflection and Transmission
    for i = 1:length(E_range)
        k0 = sqrt(E_range(i) / ekinscale);
        G = zeros(N, N);
        for m = 1:N
            for n = 1:N
                G(m, n) = exp(1i * k0 * abs(x(m) - x(n))) / (2 * 1i * k0);
            end
        end
        
        Waa = diag(U_tuned) * step / ekinscale;
        ScatMatrix = eye(N) - G * Waa;
        Phi0 = exp(1i * k0 * x);
        Phi = ScatMatrix \ transpose(Phi0);
        
        % Reflection and Transmission for tuned U1
        R_U1(j) = abs((Phi(1) - exp(1i * k0 * x(1))) / exp(1i * k0 * x(1)))^2;
        T_U1(j) = abs(Phi(end) / exp(1i * k0 * x(end)))^2;
    end
end

% Plot Transmission and Reflection for different biases
figure;
plot(U1_values, T_U1, 'r', U1_values, R_U1, 'b');
xlabel('Bias Potential [eV]');
ylabel('Coefficient');
legend('Transmission', 'Reflection');
title('Transmission and Reflection vs Bias Potential');
grid on;