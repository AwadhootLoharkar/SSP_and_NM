%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resonant Tunneling Device - Constant Background
% Num Met 4 Phys - Homework 7.8.2
% MATLAB Code
% Author: [Your Name]
% Date: [Date]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../../start-kit-student")
mystartdefaults;

% Define constants and parameters
hbar = 1.0545718e-34;   % Reduced Planck's constant (J.s)
me = 9.10938356e-31;    % Electron mass (kg)
qel = 1.60218e-19;      % Electron charge (C)
tol = 1E-12;            % Tolerance
recipunit = 1E10;       % Reciprocal unit in Å^-1
ekinscale = (hbar * recipunit)^2 / (2 * me * qel); % Conversion factor (Å^-2 to eV)

% Discretization parameters
xmin = -20; % Start of the simulation [Å]
xmax = 100; % End of the simulation [Å]
step = 0.5; % Spatial step size [Å]
x = xmin:step:xmax; % Spatial grid
N = length(x); % Number of discretization points

% Potential profile setup
U = zeros(size(x)); % Initialize potential with constant background
U(x >= 0 & x <= 15) = 0.2;    % First barrier [Å]
U(x >= 65 & x <= 80) = 0.2;   % Second barrier [Å]

% Energy parameters
E_min = 0; % Minimum energy [eV]
E_max = 0.3; % Maximum energy [eV]
dE = 0.0005; % Energy step size [eV]
E = E_min:dE:E_max; % Energy range

% Prepare matrices
k = sqrt(2 * me * E * qel) / (hbar * recipunit); % Wavevector for each energy
T = zeros(size(E)); % Transmission coefficient
R = zeros(size(E)); % Reflection coefficient
G = zeros(N, N, length(E)); % Green's function matrix for each energy

% Construct the Green's function for the reference system
for i = 1:length(E)
    k0 = sqrt(E(i) / ekinscale);
    for m = 1:N
        for n = 1:N
            G(m, n, i) = exp(1i * k0 * abs(x(m) - x(n))) / (2 * 1i * k0);
        end
    end
end

% Compute scattering matrix and solve for wavefunction
for i = 1:length(E)
    k0 = sqrt(E(i) / ekinscale);
    Waa = diag(U) * step / ekinscale; % Weighted potential matrix

    % Scattering matrix and incident wave
    ScatMatrix = eye(N) - G(:, :, i) * Waa;
    Phi0 = exp(1i * k0 * x); % Incident plane wave
    Phi = ScatMatrix \ transpose(Phi0); % Solve for scattering wavefunction
    
    % Compute transmission and reflection coefficients
    R(i) = abs((Phi(1) - exp(1i * k0 * x(1))) / exp(1i * k0 * x(1)))^2;
    T(i) = abs(Phi(end) / exp(1i * k0 * x(end)))^2;
end

% Plot Transmission and Reflection
figure;
plot(E, T, 'r', E, R, 'b');
xlabel('Energy [eV]');
ylabel('Coefficient');
legend('Transmission', 'Reflection');
title('Transmission and Reflection Coefficients for Resonant Tunneling Device');
