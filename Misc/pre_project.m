mystartdefaults

% Define the anonymous function
func = @(x) exp(-0.2*x).*cos(x);

% Create the x vector
x = 0:0.1:10;

% Evaluate the function for the x values
y = func(x);

% Generate noise
noise = rand(1, length(y)) - 0.5;

% Set the amplitude of noise
A = 1;
% Add noise to the original signal
yref = y + A * noise;

%gof(y, yref, 'short', true);

plot(x, y, 'b', x, yref, 'r')

stats = gofit(y, yref, 'tol',1e-12, 'short', false)


