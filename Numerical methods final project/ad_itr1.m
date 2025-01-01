%% Initialization
close all; clear; clc;
mystartdefaults; % Initialize constants

tic;

%% Units and Scales
recipunit = 1.0E+10; % Å^-1
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm)) / qel;

%% Grid and Potential Definition
x_step = 0.5; xmin = 0; xmax = 80; % Grid parameters (Å)
x_U = (xmin + x_step/2):x_step:(xmax - x_step/2); 
x_neg = (-20 + x_step/2):x_step:(xmin - x_step/2);
x_pos = (xmax + x_step/2):x_step:(100 - x_step/2);
xx = [x_neg, x_U, x_pos]; 

Ux = BarrierPotential(x_U, 0, 15, 0.2) + BarrierPotential(x_U, 65, 80, 0.2);
Ufull = [zeros(1, length(x_neg)), Ux, zeros(1, length(x_pos))];

E_step = 0.0005; Emin = 0.0; Emax = 0.3; 
EE = (Emin + E_step/2):E_step:(Emax - E_step/2);
damping = (hbar * 2 * pi / (1.0E-9)) / qel;

%% Biased Potential
Ebias = -0.1; 
Ubiased = Ux + Ebias * x_U / (xmax - xmin);
Ufull_biased = [zeros(1, length(x_neg)), Ubiased, Ebias * ones(1, length(x_pos))];

% Plot potential
figure;
plot(xx, Ufull, 'b', xx, Ufull_biased, 'r', 'LineWidth', 1.5);
xlabel('x [Å]'); ylabel('U(x) [eV]');
ylim([-0.15, 0.35]);
grid on;
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
legend('U(x)', 'U(x) + U_{biased}(x)', 'Location', 'Best');

%% Perturbation
Uperturbation = Ubiased - Ebias * ones(1, length(x_U));
Ucheck = [zeros(1, length(x_neg)), Uperturbation, zeros(1, length(x_pos))];

% Plot perturbation
figure;
plot(xx, Ucheck, 'b', 'LineWidth', 1.5);
xlabel('x [Å]'); ylabel('U(x) [eV]');
ylim([0, 0.35]);
grid on;

%% Wavelength and Wavevectors
commonTerm = 2 * elm * qel / hbar^2;
k0max = real(sqrt((Emax + 1i * damping) * commonTerm));
wavelength_min = 2 * pi * 1E+10 / k0max;

k1max = real(sqrt((Emax - Ebias + 1i * damping) * commonTerm));
wavelength_biased_min = 2 * pi * 1E+10 / k1max;

%% Green's Function Method
[RR, TT, AA] = RTA(Ebias, EE, damping, x_U, Uperturbation, x_step, ekinscale);

% Plot R, T, A curves
figure;
plot(-EE, RR, 'r', -EE, TT, 'g', -EE, AA, 'b', 'LineWidth', 1.5);
xlabel('E (eV)'); ylabel('R(E), T(E), A(E)');
ylim([0, 1]); grid on;
legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);

%% Tuning Bias
U0min = -0.2; U0max = 0.2;

% Plot potential for min and max biases
Umin_biased = Ux + U0min * x_U / (xmax - xmin);
Umax_biased = Ux + U0max * x_U / (xmax - xmin);

figure;
plot(xx, Ufull, 'b', xx, Umin_biased, 'r', xx, Umax_biased, 'g', 'LineWidth', 1.5);
xlabel('x [Å]'); ylabel('U(x) [eV]');
ylim([-0.25, 0.45]); grid on;
legend('U(x)', 'U(x) + Bias [U=-0.2eV]', 'U(x) + Bias [U=0.2eV]', 'Location', 'Best');

toc;
