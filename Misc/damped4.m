% Step 2: Define Initial Conditions for the Damped Vibrating String

% Given parameters
L = 0.328;              % String length [m]
M_L = 0.660E-3;         % String mass per unit length [kg/m]
F = 55.0;               % Tension force applied on string [N]
npt = 201;              % Number of spatial sample points

% Define positions of the pinches
nbPinches = 1;           % Number of pinches (e.g., one pinch in the center of the string)
pA1 = L / 2;             % First pinch position (middle of the string)
z1 = 1;                  % Amplitude scaling factor for the first pinch

% Initialize arrays for the initial conditions
f_x = zeros(npt, 1);    % Initial displacement f(x) of the string
g_x = zeros(npt, 1);    % Initial velocity g(x) of the string (set to zero initially)

% Define the x coordinate array
x = linspace(0, L, npt);  % Spatial discretization from 0 to L

% Define the initial displacement based on a single pinch in the center
for i = 1:npt
    if x(i) <= pA1
        % Linear increase in displacement up to the pinch position
        f_x(i) = (z1 * x(i)) / pA1;
    else
        % Linear decrease after the pinch position back to 0 at the string's end
        f_x(i) = (z1 * (x(i) - L)) / (pA1 - L);
    end
end

% Plot the initial displacement f(x) and initial velocity g(x)
figure;
plot(x, f_x, 'b', 'LineWidth', 2);    % Plot f(x) in blue
hold on;
plot(x, g_x, 'r--', 'LineWidth', 2);  % Plot g(x) in red dashed line
xlabel('x [m]');
ylabel('Amplitude');
title('Initial conditions f(x) and g(x)');
legend('f(x)', 'g(x)');
grid on;
hold off;

% Step 3: Compute the eigenfunctions of the vibrating string

nmax = 10; % Maximum number of modes (eigenfunctions) to compute
phi_n = zeros(npt, nmax); % Initialize array to hold eigenfunctions

% Compute eigenfunctions for n = 1, 2, ..., nmax
for n = 1:nmax
    phi_n(:, n) = sqrt(2 / L) * sin(n * pi * x / L); % Eigenfunction for mode n
end

% Plot the eigenfunctions for different n values
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

% Step 4: Calculate frequencies and periods of the eigenfunctions

% Preallocate arrays for results
omega_n = zeros(nmax, 1);  % Natural frequency (rad/s)
nu_n = zeros(nmax, 1);     % Frequency (Hz)
T_n = zeros(nmax, 1);      % Period (s)

% Calculate natural frequency (omega_n), frequency (nu_n), and period (T_n)
for n = 1:nmax
    omega_n(n) = (n * pi / L) * sqrt(F / M_L);  % Calculate natural frequency (rad/s)
    nu_n(n) = omega_n(n) / (2 * pi);            % Convert to frequency (Hz)
    T_n(n) = 1 / nu_n(n);                       % Calculate period (s)
end

% Display results in a formatted table
fprintf(' n |    nπ   |   ω_n (rad/s)  |   ν_n (Hz)   |   T_n (s)   \n');
fprintf('---|---------|----------------|--------------|--------------\n');

for n = 1:nmax
    fprintf('%2d | %8.4f | %14.4f | %13.4f | %12.4f \n', n, n * pi, omega_n(n), nu_n(n), T_n(n));
end


% Step 5: Test the orthonormality of the eigenfunctions

% Preallocate array to hold the overlap integrals
overlap_integrals = zeros(nmax, 1);

% Loop over each eigenfunction and calculate the inner product ⟨Φ_n | Φ_n⟩
for n = 1:nmax
    % Directly use the precomputed eigenfunction values for mode n
    phi_n_values = phi_n(:, n); % Extract the n-th eigenfunction
    
    % Calculate the overlap integral ⟨Φ_n | Φ_n⟩ using the trapezoidal rule
    overlap_integrals(n) = trapz(x, phi_n_values .* phi_n_values);
end

% Display the orthonormality results in a formatted table
fprintf(' n | ⟨Φ_n | Φ_n⟩ \n');
fprintf('---|---------------\n');
for n = 1:nmax
    fprintf('%2d | %13.4f \n', n, overlap_integrals(n));
end

% Step 6: Calculate overlap integrals between initial conditions and eigenfunctions

% Preallocate arrays for overlap integrals ⟨Φ_n | f⟩ and ⟨Φ_n | g⟩
f_overlap = zeros(nmax, 1);  % Overlap for f(x) (initial displacement)
g_overlap = zeros(nmax, 1);  % Overlap for g(x) (initial velocity)

% Loop to calculate overlap integrals for each mode n
for n = 1:nmax
    % Calculate the overlap integral for f(x)
    f_overlap(n) = trapz(x, phi_n(:, n) .* f_x);  % ⟨Φ_n | f⟩
    
    % Calculate the overlap integral for g(x)
    g_overlap(n) = trapz(x, phi_n(:, n) .* g_x);  % ⟨Φ_n | g⟩
end

% Display results in a formatted table
fprintf(' n | ⟨Φ_n | f⟩      | ⟨Φ_n | g⟩      \n');
fprintf('---|----------------|----------------\n');
for n = 1:nmax
    fprintf('%2d | %14.4f | %14.4f \n', n, f_overlap(n), g_overlap(n));
end

% Step 7: Compute overlap integrals ⟨Φ_n|f⟩ and ⟨Φ_n|g⟩ for the pinch condition and display table

% Preallocate arrays for the overlap integrals ⟨Φ_n|f⟩ and ⟨Φ_n|g⟩
f_overlap = zeros(nmax, 1);  % Overlap for initial displacement f(x)
g_overlap = zeros(nmax, 1);  % Overlap for initial velocity g(x)

% Loop over the modes to compute the overlap integrals using the trapezoidal rule
for n = 1:nmax
    % Compute overlap ⟨Φ_n|f⟩ for f(x) (initial displacement)
    f_overlap(n) = trapz(x, phi_n(:, n) .* f_x);
    
    % Compute overlap ⟨Φ_n|g⟩ for g(x) (initial velocity)
    g_overlap(n) = trapz(x, phi_n(:, n) .* g_x);
end

% Display the results in a table format
fprintf(' n | ⟨Φ_n | f⟩      | ⟨Φ_n | g⟩      \n');
fprintf('---|----------------|----------------\n');
for n = 1:nmax
    fprintf('%2d | %14.4f | %14.4f \n', n, f_overlap(n), g_overlap(n));
end

% Step 8: Compare f(x) and g(x) with f_approx(x) and g_approx(x)

% Compute the approximations f_approx(x) and g_approx(x)
f_approx = zeros(size(x));  % Initialize f_approx(x)
g_approx = zeros(size(x));  % Initialize g_approx(x)

% Loop over the modes to sum the contributions from each eigenfunction
for n = 1:nmax
    f_approx = f_approx + f_overlap(n) * phi_n(:, n);  % ⟨Φ_n|f⟩ * Φ_n(x)
    g_approx = g_approx + g_overlap(n) * phi_n(:, n);  % ⟨Φ_n|g⟩ * Φ_n(x)
end

% Plot comparison between f(x) and f_approx(x)
figure;
plot(x, f_x, 'b', 'LineWidth', 2); hold on;
plot(x, f_approx, 'r--', 'LineWidth', 2);
xlabel('x'); ylabel('Displacement');
legend('f(x)', 'f_{approx}(x)');
title('Comparison of f(x) and f_{approx}(x)');

% Plot comparison between g(x) and g_approx(x)
figure;
plot(x, g_x, 'b', 'LineWidth', 2); hold on;
plot(x, g_approx, 'r--', 'LineWidth', 2);
xlabel('x'); ylabel('Velocity');
legend('g(x)', 'g_{approx}(x)');
title('Comparison of g(x) and g_{approx}(x)');

% Evaluate the goodness of fit for f(x) and g(x)
gof_f = gofit(f_x, f_approx);  % Assuming gof.m is implemented
gof_g = gofit(g_x, g_approx);  % Assuming gof.m is implemented

% Display the goodness of fit
fprintf('Goodness of fit for f(x): %.4f\n', gof_f);
fprintf('Goodness of fit for g(x): %.4f\n', gof_g);

% Step 9: Compute the spatial discretizations of Ψ(x, t) for nst values of t

% Preallocate the Ψ(x, t) matrix for nst time samples
nst = 100;  % Number of time samples
tmax = 10;  % Maximum time
t = linspace(0, tmax, nst);  % Time discretization
psi_xt = zeros(length(x), nst);  % Ψ(x,t) array

% Loop over time steps to compute Ψ(x, t)
for ti = 1:nst
    % Time variable
    current_t = t(ti);
    
    % Initialize Ψ(x, t) for this time step
    psi_t = zeros(size(x));
    
    % Loop over modes to compute Ψ(x, t) as a sum over eigenfunctions
    for n = 1:nmax
        % Time evolution factor for each mode (cosine and sine terms)
        time_factor = cos(omega_n(n) * current_t) * f_overlap(n) + ...
                      (sin(omega_n(n) * current_t) / omega_n(n)) * g_overlap(n);
        
        % Add contribution of mode n to Ψ(x, t)
        psi_t = psi_t + time_factor * phi_n(:, n);
    end
    
    % Store the result for this time step
    psi_xt(:, ti) = psi_t;
end

% Plot the spatial distribution of Ψ(x, t) at different time steps
figure;
plot(x, psi_xt(:, 1), 'LineWidth', 2); hold on;
plot(x, psi_xt(:, round(nst/4)), '--', 'LineWidth', 2);
plot(x, psi_xt(:, round(nst/2)), ':', 'LineWidth', 2);
plot(x, psi_xt(:, end), '-.', 'LineWidth', 2);
xlabel('x'); ylabel('Amplitude');
legend('t=0', 't=tmax/4', 't=tmax/2', 't=tmax');
title('Spatial Distribution of Ψ(x, t) at Different Times');

% Step 10: Produce a movie of the vibrating string

% Set up video writer
video = VideoWriter('string_vibration.avi');
video.FrameRate = 10;  % Adjust frame
