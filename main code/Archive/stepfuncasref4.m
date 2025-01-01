%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step Function as Reference System with Absorption Coefficients
% Num Met 4 Phys - Homework 7.8.2 (Extended for Step Function)
% MATLAB Code
% Author: Awadhoot Loharkar
% Date: [Updated Date]
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
U0 = zeros(size(x)); % Initialize potential
x_min = 0; 
x_max = 15;

% Initial step potential
U1_initial = -0.1; % Initial step potential in eV
U0(x >= x_min & x <= x_max) = U1_initial * (1 - 0.5 * (tanh((x(x >= x_min & x <= x_max) - x_min) / 2) + 1));

% Plot the initial potential profile
figure;
plot(x, U0, 'k', 'LineWidth', 2);
xlabel('Position [Å]');
ylabel('Potential [eV]');
title('Potential Profile with Smoothed Step Function');
grid on;

%% Step 2: Define the Energy Range for Computation
E_min = 0;     % Minimum energy [eV]
E_max = 0.3;   % Maximum energy [eV]
dE = 0.5;   % Energy step size [eV]
E_range = E_min:dE:E_max;  % Energy range vector

%% Step 3: Compute Reflection, Transmission, and Absorption Coefficients
U1_values = -0.2:0.0005:0.2; % Tuning bias range
T_U1 = zeros(size(U1_values));
R_U1 = zeros(size(U1_values));
A_U1 = zeros(size(U1_values));

for j = 1:length(U1_values)
    % Update potential with current U1
    U_tuned = zeros(size(x));
    U_tuned(x >= x_min & x <= x_max) = U1_values(j) * (1 - 0.5 * (tanh((x(x >= x_min & x <= x_max) - x_min) / 2) + 1));
    
    T_avg = 0;
    R_avg = 0;
    
    % Calculate Reflection and Transmission coefficients for each energy
    for i = 1:length(E_range)
        k0 = sqrt(E_range(i) / ekinscale);
        M = eye(2);
        
        for n = 1:N
            dU = U_tuned(n) * step / ekinscale;
            M = [cos(k0 * step) - 1i * sin(k0 * step) * dU, -1i * sin(k0 * step); ...
                 -1i * sin(k0 * step) * dU, cos(k0 * step) + 1i * sin(k0 * step) * dU] * M;
        end
        
        % Reflection and Transmission coefficients
        R = abs((M(1, 2) / M(1, 1)))^2;
        T = abs(1 / M(1, 1))^2;
        
        % Accumulate coefficients
        R_avg = R_avg + R;
        T_avg = T_avg + T;
    end
    
    % Average coefficients over all energies
    R_U1(j) = R_avg / length(E_range);
    T_U1(j) = T_avg / length(E_range);
    
    % Calculate Absorption coefficient
    A_U1(j) = 1 - R_U1(j) - T_U1(j);
end

%% Step 4: Plot Reflection, Transmission, and Absorption vs Bias Potential
figure;
plot(U1_values, R_U1, 'b', 'LineWidth', 1.5); hold on;
plot(U1_values, T_U1, 'r', 'LineWidth', 1.5);
plot(U1_values, A_U1, 'g', 'LineWidth', 1.5);
xlabel('Bias Potential [eV]');
ylabel('Coefficient');
legend('Reflection R(U1)', 'Transmission T(U1)', 'Absorption A(U1)');
title('Reflection, Transmission, and Absorption vs Bias Potential');
grid on;

