%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resonant Tunneling Device Simulation with Constant Background
% Num Met 4 Phys - Homework 7.8.2
% MATLAB Code
% Author: Awadhoot Loharkar
% Date: [Date]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add necessary paths and set default configurations
addpath("../../start-kit-student");
mystartdefaults;

% Define physical constants
hbar = 1.0545718e-34;   % Reduced Planck's constant (J·s)
me = 9.10938356e-31;    % Electron mass (kg)
qel = 1.60218e-19;      % Electron charge (C)

% Numerical parameters
tol = 1E-12;            % Tolerance for numerical accuracy
recipunit = 1E10;       % Conversion factor from m^-1 to Å^-1
ekinscale = (hbar * recipunit)^2 / (2 * me * qel); % Energy scale conversion (Å^-2 to eV)

% Discretization of the spatial domain
xmin = -20; % Start position in Å
xmax = 100; % End position in Å
step = 0.05; % Spatial step size in Å
x = xmin:step:xmax; % Create spatial grid
N = length(x); % Number of grid points

% Define the potential profile (double barrier)
U = zeros(size(x)); % Initialize potential with zero background
U(x >= 0 & x <= 15) = 0.2;    % First barrier from 0 to 15 Å
U(x >= 65 & x <= 80) = 0.2;   % Second barrier from 65 to 80 Å

% Define the energy range for the incident electron
E_min = 0;      % Minimum energy in eV
E_max = 0.3;    % Maximum energy in eV
dE = 0.0005;    % Energy step size in eV
E = E_min:dE:E_max; % Generate energy array

% Preallocate arrays for transmission and reflection coefficients
T = zeros(size(E)); % Transmission coefficient
R = zeros(size(E)); % Reflection coefficient

% Precompute wavevector k for each energy
k = sqrt(2 * me * E * qel) / (hbar * recipunit); % Wavevector in Å^-1

% Green's Function initialization for reference system
G = zeros(N, N, length(E)); % Green's function matrix for each energy


%%%%%%%%%%%%%% Change the green function calculation using the eqn 7.54 %%%%%%%%%%%%%%%%%%%%%%%%%


% Construct the Green's function for the constant potential background
for i = 1:length(E)
    k0 = sqrt(E(i) / ekinscale); % Wavevector for reference system
    for m = 1:N
        for n = 1:N
            G(m, n, i) = exp(1i * k0 * abs(x(m) - x(n))) / (2 * 1i * k0);
        end
    end
end

% Loop over each energy to compute transmission and reflection
for i = 1:length(E)
    k0 = sqrt(E(i) / ekinscale); % Incident wavevector
    Waa = diag(U) * step / ekinscale; % Potential matrix scaled by step size
    
    % Construct the scattering matrix
    ScatMatrix = eye(N) - G(:, :, i) * Waa;
    
    % Define the incident plane wave (incoming wave)
    Phi0 = exp(1i * k0 * x); % Incident wavefunction
    Phi = ScatMatrix \ transpose(Phi0); % Solve for total wavefunction (scattered + incident)
    
    % Calculate reflection and transmission coefficients
    R(i) = abs((Phi(1) - exp(1i * k0 * x(1))) / exp(1i * k0 * x(1)))^2;
    T(i) = abs(Phi(end) / exp(1i * k0 * x(end)))^2;
end

% Plot the Transmission and Reflection Coefficients
figure;
plot(E, T, 'r-', 'LineWidth', 1.5); hold on;
plot(E, R, 'b--', 'LineWidth', 1.5);
xlabel('Energy [eV]', 'FontSize', 12);
ylabel('Coefficient', 'FontSize', 12);
legend('Transmission', 'Reflection', 'Location', 'best');
title('Resonant Tunneling Device - Transmission and Reflection', 'FontSize', 14);
grid on;
