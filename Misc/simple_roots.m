% Define the function f(x, c)
funct = @(x, c) exp(-x) .* (x.^2 + 2 * x + 2) - c;

% Set c = 0 and redefine func(x) as a function of x only
c = 1.1;

% Check if c is within the valid range (0, 2)
if c <= 0 || c >= 2
    disp('Root finding is not possible for c outside the range (0, 2).');
else
    % If c is valid, redefine the function as a single-variable function in x
    func = @(x) funct(x, c);
    
    % Generate x values for plotting the function f(x) over [0, 20]
    x_values = linspace(0, 20, 1000);
    y_values = arrayfun(func, x_values);
    
    % Plot the function f(x)
    %figure;
    %plot(x_values, y_values, "LineWidth", 1);
    %xlabel('x');
    %ylabel('f(x)');
    %title(['Graph of f(x) for c = ', num2str(c)]);
    %grid on;

     % Start with a crude interval [a, b]
     a = 0;
     b = 5;
     
     % Evaluate the function at the endpoints of the interval
     fa = func(a);
     fb = func(b);
     
     % Check if f(a) and f(b) have opposite signs
     if fa * fb < 0
         disp(['Crude bracketing found: f(', num2str(a), ') = ', num2str(fa), ...
               ' and f(', num2str(b), ') = ', num2str(fb)]);
         disp(['A root lies between ', num2str(a), ' and ', num2str(b)]);
     else
         disp('No sign change found in the interval [0, 20]. Trying smaller intervals...');
     end
end

% Set the desired accuracy
delta_x = 1.E-20;

% Use the bisection method to find the root within the interval [a, b]
% Assume 'bisection' is the function you already have for the bisection method
x0 = bisection( a, b,func, delta_x);

% Display the result
disp(['Approximate root x0 found: ', num2str(x0)]);
