% Vibrating String Simulation
% This script simulates the motion of a vibrating string using the wave equation
% and eigenfunction expansion method.

% 1. Define input data and discretization parameters
L = 0.328;          % String length [m]
M_L = 0.660E-3;     % String mass per unit length [kg/m]
F = 55.0;           % Tension force applied on string [N]
npt = 201;          % Number of spatial points
nst = 200;          % Number of time steps for animation (can be increased later)
nmax = 20;          % Number of eigenfunctions to use in expansion

% Calculate wave speed
c = sqrt(F / M_L);  % Wave speed [m/s]

% 2. Define the array of discretized x values
x = linspace(0, L, npt);

% 3. Compute initial conditions (pinches)
% First pinch type (single pinch at center)
f1 = zeros(1, npt);
p1 = L/2;
z1 = 0.01;  % Amplitude of the pinch
f1(x <= p1) = z1 * x(x <= p1) / p1;
f1(x > p1) = z1 * (x(x > p1) - L) / (p1 - L);

% Second pinch type (two pinches)
f2 = zeros(1, npt);
p1 = L/4; p2 = 3*L/4;
z1 = 0.01; z2 = -0.01;
f2(x <= p1) = z1 * x(x <= p1) / p1;
f2(x > p1 & x <= p2) = (z1 - z2) / (p1 - p2) * (x(x > p1 & x <= p2) - p1) + z1;
f2(x > p2) = z2 * (x(x > p2) - L) / (p2 - L);

% Initial velocity is zero for both cases
g1 = zeros(1, npt);
g2 = zeros(1, npt);

% Plot initial conditions
figure;
subplot(2,1,1);
plot(x, f1, 'b-', x, g1, 'r--');
title('Initial Conditions - Single Pinch');
legend('f(x)', 'g(x)');
xlabel('x [m]'); ylabel('Amplitude');

subplot(2,1,2);
plot(x, f2, 'b-', x, g2, 'r--');
title('Initial Conditions - Double Pinch');
legend('f(x)', 'g(x)');
xlabel('x [m]'); ylabel('Amplitude');

% 4. Compute eigenfunctions
Phi = zeros(nmax, npt);
for n = 1:nmax
    Phi(n,:) = sqrt(2/L) * sin(n*pi*x/L);
end

% Plot first few eigenfunctions
figure;
for n = 1:min(4, nmax)
    subplot(2,2,n);
    plot(x, Phi(n,:));
    title(['Eigenfunction \Phi_' num2str(n) '(x)']);
    xlabel('x [m]'); ylabel('Amplitude');
end

% 5. Display table of frequencies and periods
omega = (1:nmax)' * pi * c / L;
nu = omega / (2*pi);
T = 1 ./ nu;

disp('n    n*pi/L    omega_n    nu_n    T_n');
disp('----------------------------------------');
for n = 1:nmax
    fprintf('%2d  %8.2f  %8.2f  %6.2f  %6.4f\n', n, n*pi/L, omega(n), nu(n), T(n));
end

% 6. Test orthonormalization of eigenfunctions
overlap = zeros(nmax, nmax);
for m = 1:nmax
    for n = 1:nmax
        overlap(m,n) = trapz(x, Phi(m,:) .* Phi(n,:));
    end
end

disp('Orthonormalization test:');
disp('n    <Phi_n|Phi_n>');
disp('--------------------');
for n = 1:nmax
    fprintf('%2d  %15.12f\n', n, overlap(n,n));
end

% 7. Compute overlap integrals for initial conditions
overlap_f1 = zeros(1, nmax);
overlap_g1 = zeros(1, nmax);
overlap_f2 = zeros(1, nmax);
overlap_g2 = zeros(1, nmax);

for n = 1:nmax
    overlap_f1(n) = trapz(x, f1 .* Phi(n,:));
    overlap_g1(n) = trapz(x, g1 .* Phi(n,:));
    overlap_f2(n) = trapz(x, f2 .* Phi(n,:));
    overlap_g2(n) = trapz(x, g2 .* Phi(n,:));
end

% Display overlap integrals
disp('Overlap integrals for single pinch:');
disp('n    <Phi_n|f>    <Phi_n|g>');
disp('-------------------------------');
for n = 1:nmax
    fprintf('%2d  %10.6f  %10.6f\n', n, overlap_f1(n), overlap_g1(n));
end

disp('Overlap integrals for double pinch:');
disp('n    <Phi_n|f>    <Phi_n|g>');
disp('-------------------------------');
for n = 1:nmax
    fprintf('%2d  %10.6f  %10.6f\n', n, overlap_f2(n), overlap_g2(n));
end

% 8. Check truncation by comparing f(x) and f_approx(x)
f1_approx = zeros(1, npt);
f2_approx = zeros(1, npt);
for n = 1:nmax
    f1_approx = f1_approx + overlap_f1(n) * Phi(n,:);
    f2_approx = f2_approx + overlap_f2(n) * Phi(n,:);
end

figure;
subplot(2,1,1);
plot(x, f1, 'b-', x, f1_approx, 'r--');
title('Comparison of f(x) and f\_approx(x) - Single Pinch');
legend('f(x)', 'f\_approx(x)');
xlabel('x [m]'); ylabel('Amplitude');

subplot(2,1,2);
plot(x, f2, 'b-', x, f2_approx, 'r--');
title('Comparison of f(x) and f\_approx(x) - Double Pinch');
legend('f(x)', 'f\_approx(x)');
xlabel('x [m]'); ylabel('Amplitude');

% 9. Compute and plot string motion
t = linspace(0, 4*pi/omega(1), nst);
Psi1 = zeros(nst, npt);
Psi2 = zeros(nst, npt);

for i = 1:nst
    for n = 1:nmax
        Psi1(i,:) = Psi1(i,:) + (overlap_f1(n) * cos(omega(n)*t(i)) + overlap_g1(n)/omega(n) * sin(omega(n)*t(i))) * Phi(n,:);
        Psi2(i,:) = Psi2(i,:) + (overlap_f2(n) * cos(omega(n)*t(i)) + overlap_g2(n)/omega(n) * sin(omega(n)*t(i))) * Phi(n,:);
    end
end

% Plot string motion
figure;
subplot(2,1,1);
plot(x, Psi1);
title('String Motion - Single Pinch');
xlabel('x [m]'); ylabel('Amplitude');

subplot(2,1,2);
plot(x, Psi2);
title('String Motion - Double Pinch');
xlabel('x [m]'); ylabel('Amplitude');

% 10. Create animation
figure;
for i = 1:nst
    subplot(2,1,1);
    plot(x, Psi1(i,:), 'b-');
    title(['String Motion - Single Pinch, t = ' num2str(t(i), '%.3f') ' s']);
    xlabel('x [m]'); ylabel('Amplitude');
    axis([0 L -max(abs(f1)) max(abs(f1))]);
    
    subplot(2,1,2);
    plot(x, Psi2(i,:), 'r-');
    title(['String Motion - Double Pinch, t = ' num2str(t(i), '%.3f') ' s']);
    xlabel('x [m]'); ylabel('Amplitude');
    axis([0 L -max(abs(f2)) max(abs(f2))]);
    
    drawnow;
    pause(0.05);
end

% Note: Steps 11-13 involve changing parameters and observing results.
% These steps should be performed interactively by modifying the script
% and re-running relevant parts.

% For step 13, you can modify the initial condition f1 or f2 to be
% proportional to a specific eigenfunction, e.g.:
% m = 3;  % Choose eigenfunction
% alpha = 0.01;  % Choose amplitude
% f1 = alpha * Phi(m,:);
% Then re-run the script from step 7 onwards to observe the results.