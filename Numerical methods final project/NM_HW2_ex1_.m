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

% Plot potential
figure;
plot(xx, Ufull, 'LineWidth', 2);
xlabel('x (Å)','FontSize',26);
ylabel('U (eV)','FontSize',26);
fontsize(gca, 22,'points');
ylim([0,0.22]);
set(gca,'Box','on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
grid on;

%% Localized Green's Function Method
E_step = 0.0005; Emin = 0.0; Emax = 0.3; % Energy discretization (eV)
EE = (Emin + E_step/2):E_step:(Emax - E_step/2);
numE = length(EE);
recombT = 1.0E-9; % ns
damping = (hbar * 2 * pi / recombT) / qel;

% Calculate R(E), T(E), and A(E)
[RR, TT, AA] = RTA(0, EE, damping, x_U, Ux, x_step, ekinscale);

% Plot A(E)
figure;
plot(EE, AA, 'LineWidth', 2);
xlabel('E (eV)','FontSize',26);
ylabel('A(E)','FontSize',26);
fontsize(gca, 22,'points');
ylim([0,0.03]);
set(gca,'Box','on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
grid on;

%% Resonant Tunneling Diode (RTD) Model
resonanceE = [0.01075, 0.04325, 0.09525, 0.16325, 0.23975];
nonresonantE = (resonanceE(3) + resonanceE(4)) / 2;
resoE = [resonanceE, nonresonantE];
[~, numResonantEnergies] = size(resoE);

wavefunctions = zeros(numResonantEnergies, length(xx));
potentialMatrix = zeros(length(x_U), length(x_U));

for energyIndex = 1:numResonantEnergies
    waveNumber = sqrt((resoE(energyIndex) + 1i * damping) / ekinscale);
    
    for potentialIndex = 1:length(x_U)
        waveFunctionAtWaveNumber(potentialIndex) = exp(1i * waveNumber * x_U(potentialIndex));
        potentialMatrix(potentialIndex, potentialIndex) = x_step * Ux(potentialIndex) / ekinscale;
        
        for greenIndex = 1:length(x_U)
            greenMatrix(potentialIndex, greenIndex) = GreensFun(0, x_U(potentialIndex), x_U(greenIndex), resoE(energyIndex), damping, ekinscale);
        end
    end
    
    scatteringMatrix = eye(length(x_U)) - greenMatrix * potentialMatrix;
    waveFunctionSolution = scatteringMatrix \ (waveFunctionAtWaveNumber.'); 
    
    for positionIndex = 1:length(xx)
        extraTermValue = ExtraTerm(xx(positionIndex), x_U, Ux, x_step, waveFunctionSolution, 0, resoE(energyIndex), damping, ekinscale);
        
        if energyIndex < numResonantEnergies
            wavefunctions(energyIndex, positionIndex) = exp((xx(positionIndex)) * 1i * waveNumber) + extraTermValue;
        else
            nonResonantWaveFunction(positionIndex) = exp((xx(positionIndex)) * 1i * waveNumber) + extraTermValue;
        end
    end
end

resonantProbabilities = abs(wavefunctions).^2;
nonResonantProbability = abs(nonResonantWaveFunction).^2;

% Plot scaled resonant wavefunctions
figure;
hold on;
plot(xx, (100/0.2)*Ufull, 'LineWidth', 2, 'Color', 'black');
plot(xx, (100/4.)*resonantProbabilities(1,:), 'LineWidth', 1);
plot(xx, (100/27.)*resonantProbabilities(2,:), 'LineWidth', 1);
plot(xx, (100/112.)*resonantProbabilities(3,:), 'LineWidth', 1);
plot(xx, (100/30.)*resonantProbabilities(4,:), 'LineWidth', 1);
plot(xx, (100/6.)*resonantProbabilities(5,:), 'LineWidth', 1);
plot(xx, (100/4.)* nonResonantProbability, 'LineWidth', 2, 'Color', [0 0 0]+0.5);
xlabel('x (Å)','FontSize',26);
ylabel('|\psi(x)|^2','FontSize',26);
fontsize(gca, 22,'points');
legend('U(x)', '|\psi_1(x)|^2', '|\psi_2(x)|^2', '|\psi_3(x)|^2', '|\psi_4(x)|^2', '|\psi_5(x)|^2', '|\psi_{nr}(x)|^2', 'Location', 'eastoutside');
ylim([0,125]);
set(gca,'Box','on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
grid on;

%% Crude Approach to an Applied Bias
Ux = BarrierPotential(x_U, 0, 15, 0.2) + BarrierPotential(x_U, 65, 80, 0.1);
[RR, TT, AA] = RTA(0, EE, damping, x_U, Ux, x_step, ekinscale);

% Plot A(E) for modified quantum junction
figure;
plot(EE, AA, 'LineWidth', 2);
xlabel('E (eV)','FontSize',26);
ylabel('A(E)','FontSize',26);
fontsize(gca, 22,'points');
ylim([0,0.005]);
set(gca,'Box','on', 'LineWidth', 1, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
set(gcf, 'Color', [0.7 0.7 0.7]);
grid on;

toc;