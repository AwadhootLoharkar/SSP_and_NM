% Given parameters
L = 0.328;              % String length [m]
M_L = 0.660E-3;         % String mass per unit length [kg/m]
F = 55.0;               % Tension force applied on string [N]
npt = 201;              % Number of spatial sample points
nst = 200;               % Number of time sample points
damp_factor = 200;      % Damping factor

% Define positions of the pinches
nbPinches = 1;         % Number of pinches
pA1 = L / 2;           % First pinch position
pB1 =  L / 4;          % Second pinch position
pB2 = 3 * L / 4;       % Second pinch position
z1 = 1;                % Amplitude scaling factor for the first pinch
z2 = -1;               % Amplitude scaling factor for the second pinch (opposite direction)

% Initialize arrays for the initial conditions
f_x = zeros(npt, 1);   % Initial displacement
g_x = zeros(npt, 1);   % Initial velocity

% Define the x coordinate array
x = linspace(0, L, npt);  % Spatial discretization from 0 to L

% Define the initial displacement based on the pinch functions
if nbPinches == 1
    for i = 1:npt
        % First pinch at L/4 (positive direction)
        if x(i) <= pA1
            f_x(i) = (z1 * x(i)) / pA1;  % Pinch increases linearly to p1
        else 
            f_x(i) = (z1 * (x(i) - L)) / (pA1 - L); % Decrease towards 0 until p2
        end
    end    
else
    for i = 1:npt
        if x(i) <= pB1
            f_x(i) = (z1 * x(i)) / pB1;
        elseif x(i) <= pB2
            f_x(i) = ((z1 - z2)*(x(i) - pB1)/ (pB1 - pB2)) + z1;
        elseif x(i) > pB2
            f_x(i) = (z2 * (x(i) - L)) / (pB2 - L);
        end
    end
end

% Compute and plot f(x) and g(x)
figure;
plot(x, f_x, 'b', 'LineWidth', 2);
hold on;
plot(x, g_x, 'r--', 'LineWidth', 2);
xlabel('x [m]');
ylabel('Amplitude');
title('Initial conditions f(x) and g(x) with two opposite pinches');
legend('f(x)', 'g(x)');
grid on;
hold off;

% Step 4: Compute the numerical arrays describing the eigenfunctions.
% Compute eigenfunctions for n = 1, 2, ..., nmax
nmax = 10; % Define the maximum number of modes to calculate
phi_n = zeros(npt, nmax); % Initialize an array to hold eigenfunctions

for n = 1:nmax
    phi_n(:, n) = sqrt(2/L) * sin(n * pi * x / L); % Compute eigenfunction
end

% Plot the eigenfunctions for various n values
figure;
hold on;
for n = 1:nmax
    plot(x, phi_n(:, n), 'DisplayName', sprintf('\\Phi_{%d}(x)', n), 'LineWidth', 1.5);
end
xlabel('x [m]');
ylabel('Eigenfunctions \Phi_n(x)');
title('Eigenfunctions of the Vibrating String');
legend show; % Display the legend for each mode
grid on;
hold off;

% Step 5: Calculate the overlap integrals
% Calculate natural frequency (omega_n), frequency (nu_n), and period (T_n)
omega_n = zeros(nmax, 1);  % Preallocate omega_n array
for n = 1:nmax
    omega_n(n) = (n * pi / L) * sqrt(F / M_L);  % Calculate omega_n
    nu_n(n) = omega_n(n) / (2 * pi);             % Calculate frequency nu_n
    T_n(n) = 1 / nu_n(n);
end

% Display results in a formatted table
fprintf(' n |    nπ   |   ω_n (rad/s)  |   ν_n (Hz)   |   T_n (s)   \n');
fprintf('---|---------|----------------|--------------|--------------\n');

for n = 1:nmax
    fprintf('%2d | %8.4f | %14.4f | %13.4f | %12.4f \n', ...
        n, n_pi(n), omega_n(n), nu_n(n), T_n(n));
end

% Step 6: Test the orthonormalization of eigenfunctions

% Preallocate an array to hold the overlap integrals
overlap_integrals = zeros(nmax, 1);

% Loop to calculate the overlap integral for each n
for n = 1:nmax
    % Directly access the precomputed eigenfunction values
    phi_n_values = phi_n(:, n);  % Use the n-th eigenfunction directly
    
    % Calculate the overlap integral ⟨Φ_n | Φ_n⟩ using the trapezoidal rule
    overlap_integrals(n) = trapz(x, phi_n_values .* phi_n_values);
end

% Display results in a formatted table
fprintf(' n | ⟨Φ_n | Φ_n⟩ \n');
fprintf('---|---------------\n');

for n = 1:nmax
    fprintf('%2d | %13.4f \n', n, overlap_integrals(n));
end

% Calculate overlap integrals between initial conditions and eigenfunctions
f_overlap = zeros(nmax, 1);  % ⟨Φ_n | f⟩
g_overlap = zeros(nmax, 1);  % ⟨Φ_n | g⟩
for n = 1:nmax
    f_overlap(n) = trapz(x, phi_n(:, n) .* f_x);  % ⟨Φ_n | f⟩
    g_overlap(n) = trapz(x, phi_n(:, n) .* g_x);  % ⟨Φ_n | g⟩
end
% Display results in a formatted table
fprintf(' n | ⟨Φ_n | f⟩      | ⟨Φ_n | g⟩      \n');
fprintf('---|----------------|----------------\n');

for n = 1:nmax
    fprintf('%2d | %14.4f | %14.4f \n', n, f_overlap(n), g_overlap(n));
end

% Define time discretization over one period with damping
t = linspace(0, 2 * pi / omega_n(1), nst); % Time array for one period with nst points
gamma = damp_factor; % Damping factor
psi_xt = zeros(npt, nst); % Initialize the array for Ψ(x, t)

% Compute Ψ(x, t) with damping for each time sample
for j = 1:nst
    for n = 1:nmax
        % Exponential decay factor e^(-gamma * t(j))
        damping = exp(-gamma * t(j));
        
        % Compute the contribution of each eigenfunction to Ψ(x, t) with damping
        psi_xt(:, j) = psi_xt(:, j) + damping * ...
            (f_overlap(n) * cos(omega_n(n) * t(j)) + g_overlap(n) * sin(omega_n(n) * t(j))) * phi_n(:, n);
    end
end

% Plot the amplitude patterns at specific time instances with damping
figure;
hold on;
for j = 1:5:nst % Plot every 5th time step for clarity
    plot(x, psi_xt(:, j), 'DisplayName', sprintf('t = %.2f s', t(j)), 'LineWidth', 1.5);
end
xlabel('x [m]');
ylabel('Ψ(x, t)');
title('Amplitude Patterns of the Damped Vibrating String at Different Times');
legend show;
grid on;
hold off;

% Produce a movie of the string vibrating with damping during at least 4π/ω_1
n_frames = 200;  % Number of frames for the movie
time_duration = 4 * pi / omega_n(1); % Duration to observe the motion
t_movie = linspace(0, time_duration, n_frames); % Time array for the movie
frames(n_frames) = struct('cdata', [], 'colormap', []);

figure;
axis([0 L -1 1]); % Set axis limits
xlabel('x [m]');
ylabel('Ψ(x, t)');
title('Damped Vibrating String Motion');
grid on;
hold on;

for j = 1:n_frames
    cla; % Clear previous plot
    
    current_psi = zeros(npt, 1);
    for n = 1:nmax
        damping = exp(-gamma * t_movie(j));
        current_psi = current_psi + damping * ...
            (f_overlap(n) * cos(omega_n(n) * t_movie(j)) + g_overlap(n) * sin(omega_n(n) * t_movie(j))) * phi_n(:, n);
    end
    
    plot(x, current_psi, 'LineWidth', 2);
    ylim([-1 1]); % Keep consistent plot limits
    title(sprintf('Damped String at t = %.2f s', t_movie(j)));
    
    frames(j) = getframe(gcf); % Capture frame
end

% Save the movie with damping
movie_name = 'damped_vibrating_string_movie.mp4';
v = VideoWriter(movie_name, 'MPEG-4');
open(v);
writeVideo(v, frames);
close(v);
