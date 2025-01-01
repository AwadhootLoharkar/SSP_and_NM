%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%    Num Met 4 Phys - Ex 7.8.1  %
%          15,22 Oct 24         %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath("../../start-kit-student")

mystartdefaults;
tol = 1E-12;
recipunit = 1E10; % Length of the Brillouin zone [m^-1]
ekinscale = (hbar*recipunit)^2/(2*elm*qel); % Convertion factor between square of wavevector [A^-2] to energy [eV]

%% Input parameters

% Points for discretization
xminA = 0; % Start of the perturbation [A]
xmaxA = 5; % End of the perturbation [A]

xminB = -5;% Start of the simulation [A]
xmaxB = 10; % End of the simulation [A]

step = 0.25; % Same step everywhere [A]

dimA = (xmaxA - xminA)/step;
dimB = (xmaxB - xminB)/step;

% Energies
Uzero = 2; % Height of the step [eV]
Ezero = 1; % Energy of the incident electron [eV]


%% Points for discretization

for i=1:dimA
    xA(i) = xminA+step/2+step*(i-1);
end

for i=1:dimB
    xB(i) = xminB+step/2+step*(i-1);
end


%% Checking min and max energies

Emax = ekinscale*(xmaxA-xminA);
Emin = ekinscale*step;

fprintf("INFO: Minimum energy allowed is %d eV\n", Emin);
fprintf("INFO: Maximum energy allowed is %d eV\n", Emax);

if(Uzero > Emax || Uzero < Emin)
    fprintf("WARN: Energy of the step will be affected by aliasing.\n");
end
if(Ezero > Emax || Ezero < Emin)
    fprintf("WARN: Energy of the incident electron will be affected by aliasing.\n");
end

%% Spatial part of the incident plane wave inside the perturbation

kzero = sqrt(Ezero/ekinscale);
PhiAzero = exp(1i.*kzero.*xA); % eq 7.15


%% Perturbation

Va = ones(1, dimA).*Uzero; % [eV]
Waa = diag(Va)*step/ekinscale; % eq 7.9 [A^-2]


%% Plot of the perturbation

Vb = zeros(1, dimB); j=1;
for i=1:dimB
    if(xB(i) >= xminA && xB(i) <= xmaxA)
        Vb(i) = Va(j); j = j+1; % Just for plots
    end
end

figV = figure('name', 'Perturbation', 'NumberTitle', 'off');
hold on;
plot(xB, Vb);
xlabel("x [A]");
ylabel("V [eV]");
movegui(figV, 'northwest');


%% Green's function

Gaa = zeros(dimA, dimA);
for i=1:dimA
    for j=1:dimA
        Gaa(i, j) = exp(1i*kzero*abs(xA(i)-xA(j)))/(2*1i*kzero); % eq 7.34 [A]
    end
end


%% Scattering matrix

ScatMatrix = eye(dimA, dimA)-Gaa*Waa; % Scattering matrix eq 7.56
PhiA = ScatMatrix\transpose(PhiAzero); % eq 7.56


%% Squared modulus of the wavefunction inside the perturbation

PhiASquared = abs(PhiA).^2;


%% Plot of the wavefunction inside the perturbation

figPhiA = figure('name', 'Inside the perturbation', 'NumberTitle', 'off');
hold on;
plot(xA, PhiASquared);
xlabel("x [A]");
ylabel("|\phi_A|^2");
movegui(figPhiA, 'north');


%% Compute the full wavefunction

PhiB = exp(1i*kzero.*xB); % PhiBZero
for b=1:dimB % eq 7.54 (Huygens' principle)
    for a=1:dimA
        PhiB(b) = PhiB(b) + exp(1i*kzero*abs(xB(b)-xA(a)))/(2*1i*kzero)*Waa(a, a)*PhiA(a);
    end
end


%% Squared modulus of the full wavefunction

PhiBSquared = abs(PhiB).^2;

%% Plot of the full wavefunction

figPhiB = figure('name', 'Full simulation', 'NumberTitle', 'off');
hold on;
plot(xB, PhiBSquared);
xlabel("x [A]");
ylabel("|\phi_B|^2");
movegui(figPhiB, 'northeast');

%% Reflection & transmission coefficients

R = abs((PhiB(1)-exp(1i*kzero*xB(1)))/exp(1i*kzero*xB(1)))^2; % PhiBScattered/PhiBZero in an arbitrary point before
T = abs(PhiB(end)/exp(1i*kzero*xB(end)))^2; % PhiB/PhiAZero in an aribtrary point after

fprintf("\nRESULT: Reflection coefficient: %d\n", R);
fprintf("RESULT: Transmission coefficient: %d\n", T);

if(abs(R+T-1) > tol)
    fprintf("WARN: The sum of reflection and transmission coefficients is not 1: %d\n", R+T);
end


%% Probability current

speedFactor = recipunit*hbar/elm; % Factor to get the current probability in m/s
speedFactor = speedFactor*10e-6; % Factor to get the current probability in M m/s

J = zeros(1, dimB);
for i=1:dimB-1
    grad = (PhiB(i+1)-PhiB(i))/step; % 2 points approximation
    J(i) = speedFactor*real(-1i*conj(PhiB(i))*grad); % eq 7.91
end
J(end) = J(end-1);


%% Plot of the probability current

figJ = figure('name', 'Probability current', 'NumberTitle', 'off');
hold on;
plot(xB, J);
xlabel("x [A]");
ylabel("J [10^6 m/s]");
ylim([0 2*max(J)]);
movegui(figJ, 'southwest');

%% Question 11 skipped
