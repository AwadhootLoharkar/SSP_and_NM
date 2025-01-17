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
numE = length(EE);
recombT = 1.0E-9; % ns
damping = (hbar * 2 * pi / recombT) / qel;

%% Biased Potential
Ebias = -0.1; 
Ubiased = Ux + Ebias * x_U / (xmax - xmin);
Ufull_biased = [zeros(1, length(x_neg)), Ubiased, Ebias * ones(1, length(x_pos))];

% Plot potential
figure;
plot(xx, Ufull, 'b', xx, Ufull_biased, 'r');
xlabel('x [Å]');
ylabel('U(x) [eV]');
ylim([-0.15, 0.35]);
fontsize(gca, 22, 'points');
set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
legend('U(x)', 'U(x) + U_{biased}(x)', 'Location', 'Best');
grid on;

%% Perturbation
Uperturbation = Ubiased - Ebias * ones(1, length(x_U));
Ucheck = [zeros(1, length(x_neg)), Uperturbation, zeros(1, length(x_pos))];

% Plot perturbation
figure;
plot(xx, Ucheck, 'b');
xlabel('x [Å]');
ylabel('U(x) [eV]');
ylim([0, 0.35]);
fontsize(gca, 22, 'points');
set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
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
plot(-EE, RR, 'r', -EE, TT, 'g', -EE, AA, 'b', 'LineWidth', 2);
xlabel('E (eV)');
ylabel('R(E), T(E), A(E)');
ylim([0, 1]);
fontsize(gca, 22, 'points');
set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');
grid on;

%% Tuning Bias
U0min = -0.2; U0max = 0.2;

Umin_biased = Ux + U0min * x_U / (xmax - xmin);
Umax_biased = Ux + U0max * x_U / (xmax - xmin);
U0min_biased = [zeros(1, length(x_neg)), Umin_biased, U0min * ones(1, length(x_pos))];
U0max_biased = [zeros(1, length(x_neg)), Umax_biased, U0max * ones(1, length(x_pos))];

% Plot potential for min and max biases
figure;
plot(xx, Ufull, 'b', xx, U0min_biased, 'r', xx, U0max_biased, 'g');
xlabel('x [Å]');
ylabel('U(x) [eV]');
ylim([-0.25, 0.45]);
fontsize(gca, 22, 'points');
set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
legend('U(x)', 'U(x) + Bias [U=-0.2eV]', 'U(x) + Bias [U=0.2eV]', 'Location', 'Best');
grid on;

%% RTA for the Electron at the Fermi Level
setE = 0.01; % eV, E0
EE = (Emin + E_step/2 + U0min):E_step:(Emax - E_step/2 + U0min); % eV, Discretized energy

for i = 1:length(EE)
    for j = 1:length(x_U)
        perturbation(j) = EE(i) * x_U(j) / (xmax - xmin) - EE(i);
    end
    Upertu = Ux + perturbation;
    [RR(i), TT(i), AA(i)] = RTA_iter(EE(i), setE, damping, x_U, Upertu, x_step, ekinscale);
end

% Plot R(E), T(E), A(E)
figure;
plot(-EE, RR, 'LineWidth', 2);
hold on;
plot(-EE, TT, 'LineWidth', 2);
plot(-EE, AA, 'LineWidth', 2);
xlabel('E (eV)');
ylabel('R(E), T(E), A(E)');
fontsize(gca, 22, 'points');
ylim([0, 1]);
set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');
grid on;

% Plot R(E), T(E), A(E) for a specific range
figure;
plot(-EE, RR, 'LineWidth', 2);
hold on;
plot(-EE, TT, 'LineWidth', 2);
plot(-EE, AA, 'LineWidth', 2);
xlabel('E (eV)');
ylabel('R(E), T(E), A(E)');
fontsize(gca, 22, 'points');
xlim([0.06, 0.07]);
ylim([0, 0.2]);
set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');
grid on;

toc;