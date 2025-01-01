%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%   Num Met 4 Phys - Ex 2.13.6  %
%           24 Sept 24          %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../../start-kit-student")
addpath("../bin")

mystartdefaults;
tol = 1E-12;

%% Input data

L = 0.328; % Length of the string [m]
rho = 0.660E-3; % String mass par unit length [kg/m]
F = 55.0; % Tension force applied on string [N]
C = sqrt(F/rho); % Speed of the wave

npt = 201; % Number of sample points for space
nmax = 30; % Number of terms of the (in)finite sum
nst = 10; % Step increment for time
tmax = 1000; % End of the time simulation (start is 0)

if (nmax >= npt)
    warning("The number of modes is greater than the number of sample points for space.");
end


%% Defining useful arrays

x = linspace(0, L, npt); % Space discretization of the string
f1 = OnePinch(x, x(end)/2, 2);
f2 = TwoPinches(x, x(end)/4, 3*x(end)/4, 2, -2);
%f2 = TwoPinches(x, x(end)/7, 3*x(end)/4, 2, -2);

f = f2;
g = x.*0;

% Plots

fig1 = figure('NumberTitle', 'off', 'name','Useful arrays: Pinches'); hold on;
plot(x, f1);
plot(x, f2);
legend("One Pinch", "Two Pinches");
xlim([0, x(end)]);
movegui(fig1, 'northwest');

%% Eigenfunctions

Phi = zeros(nmax, npt);
for n=1:nmax
    Phi(n, :) = sqrt(2/L).*sin(n*pi*x/L); % eq 2.48
end

% Plots

fig2 = figure('NumberTitle', 'off', 'name','Eingenfunctions'); hold on;
plot(x, Phi);
xlim([0, x(end)]);
movegui(fig2, 'north');


%% Table

Table = zeros(nmax, 5);
N = 1:nmax;               % n
waveVector = N.*pi/L;     % k = n*pi/L
omega = waveVector.*C;        % omega_n = n*pi*c/L
nu = omega./(2*pi);   % nu_n = omega_n/2*pi
T = 1./nu;        % T_n = 1/nu_n

printtable([N ; waveVector ; omega ; nu ; T]', 'ColName', {'n', 'n*pi/L', 'omega_n', 'nu_n', 'T_n'}, 'Normalized', true);

% Fundamental frequency is 440Hz, it's an A.


%% Orthonormalization

fprintf("***** Orthonormalization check *****\n\n")

Ortho = zeros(nmax, nmax);

for i=1:nmax
    for j=1:nmax
        Ortho(i, j) = trapz(x, conj(Phi(i, :)).*Phi(j, :)); % eq 2.51
    end
end
Ortho(abs(Ortho)<tol) = 0;

printtable(Ortho);


fprintf("***** Overlaps *****\n\n")

OverlapF = zeros(1, nmax);
OverlapG = zeros(1, nmax);

for i=1:nmax
    OverlapF(i) = trapz(x, conj(Phi(i, :)).*f); % eq 2.57
    OverlapG(i) = trapz(x, conj(Phi(i, :)).*g); % eq 2.57
end
OverlapF(abs(OverlapF)<tol) = 0;
OverlapG(abs(OverlapG)<tol) = 0;

printtable([(1:n) ; OverlapF ; OverlapG]', 'ColName', {'n', '<phi_n|f>', '<phi_n|g>'}, 'Normalized', true, 'Integer', 1);


%% Computing approximated initial functions f & g

fapprox = zeros(1, npt);

for i=1:length(x)
    fapprox(i) = sum(OverlapF'.*Phi(:, i)); % eq 2.65
end

gapprox = zeros(1, npt);

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

%% Some frames

fig4 = figure('name', 'Some frames - Damped string', 'NumberTitle', 'off');
hold on;

time = linspace(0, max(T), 20);

for t=time
    Psi = zeros(1, length(x));
    for n=1:nmax
        Psi = Psi + Phi(n, :).*(OverlapF(n).*cos(omega(n)*t)+1/omega(n)*OverlapG(n)*sin(omega(n)*t));
    end
    plot(x, Psi);
end

xlim([0, x(end)]);
ylim([-max(f)-1, max(f1)+1]);
movegui(fig4, 'southeast');

%% Movie

Psi = zeros(1, length(x));

fig5 = figure('name', 'Movie - Damped string', 'NumberTitle', 'off');
hold on;
plot(x, Psi, 'XDataSource', 'x', 'YDataSource', 'Psi');
xlim([0, x(end)]);
ylim([-max(f)-1, max(f1)+1]);
movegui(fig5, 'southwest');

for t=0:nst:tmax
    Psi = zeros(1, length(x));
    for n=1:nmax
        Psi = Psi + Phi(n, :).*(OverlapF(n).*cos(omega(n)*t)+1/omega(n)*OverlapG(n)*sin(omega(n)*t));
    end
    refreshdata();
    drawnow();
end

