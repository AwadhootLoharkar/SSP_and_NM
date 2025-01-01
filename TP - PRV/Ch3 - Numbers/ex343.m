%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%    Num Met 4 Phys - Ex 3.4.3  %
%           11 Oct 24           %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../../start-kit-student")
addpath("../bin")

mystartdefaults;
tol = 1E-12;

%% Question 1+2

c = 10^666;
fprintf("INFO: Findind the root of f(x)=exp(-x)*(x^2+2*x+2)-%24.16E\n", c);

funct = @(x, c) exp(-x).*(x.^2+2.*x+2)-c; % Two variables function
func = @(x) funct(x, c);

x = 0:0.1:20;

fig1 = figure('name', 'Plots of the function', 'NumberTitle', 'off');
hold on; ylim([-2, 2]);

for c_temp = [-1, 0, 1.E-6, 0.1 0.5, 1, 1.99, 2, 2.01, 3]
    y = funct(x, c_temp);
    plot(x, y);
end
legend("-1", "0", "1.E-6", "0.1", "0.5", "1", "1.99", "2", "2.01", "3");
plotzeros(); xlabel("x"); ylabel("f(x)");

%% Question 3

if(c >= 2 || c < 0)
    error("ERROR: The function does not have any root.");
end

%% Question 4

b = realmin;
while(func(b) >= 0)
    b = 2*b;
end
a = b/2;
fprintf("INFO: Bracketing in interval [%24.16E, %24.16E]\n", a, b);

%% Question 5

dx = 1.E-6;
[approx, iter] = bisection (a, b, func, dx);
fprintf("INFO: Result given with precision %24.16E (it needed %d iterations)\n", dx, iter);
fprintf("The root is approximatly x_0=%24.16E\n", approx);
fprintf("f(x_0)=%24.16E\n", func(approx));

%% Question 6

fprintf("\n");

a = approx-dx;
b = approx+dx;
fprintf("INFO: Bracketing in interval [%24.16E, %24.16E]\n", a, b);
dx = eps(approx); % Best accuracy around approx

[approx, iter] = bisection (a, b, func, dx);
fprintf("INFO: Result given with precision %24.16E (it needed %d iterations)\n", dx, iter);
fprintf("The root is approximatly x_0=%24.16E\n", approx);
fprintf("f(x_0)=%24.16E\n", func(approx));

%% Question 7

% It's not looping forever. The algorithm stops when the length of the
% bracketing interval is not greater than dx. At some point, the length
% will become zero (after being a subnormal number...)

%% Question 8

