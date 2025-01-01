%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%   Num Met 4 Phys - Ex 2.14.5  %
%           1 Oct 24            %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../../start-kit-student")
addpath("../bin")

mystartdefaults;
tol = 1E-12;

%% Input data (Q1)

L = 0.328; % Length of the string [m]
rho = 0.660E-3; % String mass par unit length [kg/m]
F = 55.0; % Tension force applied on string [N]
C = sqrt(F/rho); % Speed of the wave
gamma = 100; % Damping factor [Hz]
resolution = 20; % Frames per period

npt = 501; % Number of sample points for space (minus 1, due to the discretization definition)
nmax = 10; % Number of terms of the (in)finite sum
tmax = 2*log(10)/gamma; % End of the time simulation to observe a 100x decay

if (nmax >= npt)
    warning("The number of modes is greater than the number of sample points for space.");
end


%% Defining useful arrays (Q2, Q3)


xstep = L/npt; % In A
for i=1:npt-1
    x(i) = xstep/2+xstep*(i-1); % Space discretization of the string
end

f1 = OnePinch(x, x(end)/2, 2);
f2 = TwoPinches(x, x(end)/4, 3*x(end)/4, 2, -2);
%f2 = TwoPinches(x, x(end)/7, 3*x(end)/4, 2, -2);

f = f2; % Initial condition (position) used in the following
g = x.*0; % Initial condition (speed) used in the following

% Plots

fig1 = figure('NumberTitle', 'off', 'name','Useful arrays: Pinches'); hold on;
plot(x, f1);
plot(x, f2);
legend("One Pinch", "Two Pinches");
xlim([0, x(end)]);
movegui(fig1, 'northwest');

%% Eigenfunctions (Q4)

Phi = zeros(nmax, npt-1);
for n=1:nmax
    Phi(n, :) = sqrt(2/L).*sin(n*pi*x/L); % eq 2.48
end

% Plots

fig2 = figure('NumberTitle', 'off', 'name','Eingenfunctions'); hold on;
plot(x, Phi);
xlim([0, x(end)]);
movegui(fig2, 'north');


%% Table (Q5)

Table = zeros(nmax, 5);
N = 1:nmax;               % n
waveVector = N.*pi/L;     % k = n*pi/L
omega = waveVector.*C;        % omega_n = n*pi*c/L
nu = omega./(2*pi);   % nu_n = omega_n/2*pi
T = 1./nu;        % T_n = 1/nu_n

printtable([N ; waveVector ; omega ; nu ; T]', 'Integer', [1], 'Precision', [3], 'ColName', {'n', 'n*pi/L', 'omega_n', 'nu_n', 'T_n'}, 'Normalized', true);

% Fundamental frequency is 440Hz, it's an A.

%% Orthonormalization (Q6, Q7)

fprintf("***** Orthonormalization check *****\n\n")

Ortho = zeros(nmax, nmax);

for i=1:nmax
    for j=1:nmax
        Ortho(i, j) = trapz(x, conj(Phi(i, :)).*Phi(j, :)); % eq 2.51
    end
end
Ortho(abs(Ortho)<tol) = 0;

printtable(Ortho, 'Precision', [3]);


fprintf("***** Overlaps *****\n\n")

OverlapF = zeros(1, nmax);
OverlapG = zeros(1, nmax);

for i=1:nmax
    OverlapF(i) = trapz(x, conj(Phi(i, :)).*f); % eq 2.57
    OverlapG(i) = trapz(x, conj(Phi(i, :)).*g); % eq 2.57
end
OverlapF(abs(OverlapF)<tol) = 0;
OverlapG(abs(OverlapG)<tol) = 0;

printtable([(1:n) ; OverlapF ; OverlapG]', 'Precision', [3], 'ColName', {'n', '<phi_n|f>', '<phi_n|g>'}, 'Normalized', true, 'Integer', 1);


%% Computing approximated initial functions f & g (Q8)

fapprox = zeros(1, npt-1);

for i=1:length(x)
    fapprox(i) = sum(OverlapF'.*Phi(:, i)); % eq 2.65
end

gapprox = zeros(1, npt-1);

for i=1:length(x)
    gapprox(i) = sum(OverlapG'.*Phi(:, i)); % eq 2.65
end

% Plot

fig3 = figure('name', 'Approximated initial functions', 'NumberTitle', 'off');
hold on;
plot(x, fapprox);
plot(x, gapprox);
legend("f_{approx}", 'g_{approx}');
xlim([0, x(end)]);
movegui(fig3, 'northeast');

gof(f, fapprox)

%% Some frames (Q9)

fig4 = figure('name', 'Some frames - Undamped string', 'NumberTitle', 'off');
hold on;

time = linspace(0, T(1), resolution);

for t=time
    Psi = zeros(1, length(x));
    for n=1:nmax
        Psi = Psi + Phi(n, :).*(OverlapF(n).*cos(omega(n)*t)+1/omega(n)*OverlapG(n)*sin(omega(n)*t)).*exp(-1*gamma*t);
    end
    plot(x, Psi);
end

xlim([0, x(end)]);
ylim([-max(f)-1, max(f1)+1]);
movegui(fig4, 'southeast');

%% Movie (Q10)

Psi = zeros(1, length(x)); % Contains 1 wavefunction at a specific time, then replaced by the next one
nst = T(1)/resolution; % Step increment for time

fig5 = figure('name', 'Movie - Undamped string', 'NumberTitle', 'off');
hold on;
plot(x, Psi, 'XDataSource', 'x', 'YDataSource', 'Psi');
xlim([0, x(end)]);
ylim([-max(f)-1, max(f1)+1]);
movegui(fig5, 'southwest');

for t=0:nst:tmax
    Psi = zeros(1, length(x));
    for n=1:nmax
        Psi = Psi + Phi(n, :).*(OverlapF(n).*cos(omega(n)*t)+1/omega(n)*OverlapG(n)*sin(omega(n)*t)).*exp(-1*gamma*t);
    end
    refreshdata();
    drawnow();
end