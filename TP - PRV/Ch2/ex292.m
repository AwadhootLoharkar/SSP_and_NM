% Euler.m
clear all;   % Clear all variables/functions in memory
close all;   % Close any previous figure
clc;         % Clear terminal
addpath("../../start-kit-student")

pi           % pi is predefined in Matlab/Octave
1i           % Imaginary number in Matlab/Octave
graph=true;  % Logical variable

fmt=2;       % Values to test = 1,2,3
switch fmt   % Testing numerical variables 
  case 1
    fmt_string = '%d %d %d \n';
  case 2
    fmt_string = '%15.6E %15.6E %15.6E \n'; 
  case 3
    fmt_string = '%#15.6G %#15.6G %#15.6G \n';
  otherwise
    error('fmt value not recognized')
end
    
n = 101 ;   xmin = 0; xmax = 2*pi; step = (xmax-xmin)/(n-1);
for k=1:n
  x(k) = xmin + (k-1) * step;  % Compute abscissas
end

z = exp(1i*x); % Compute ordinates

f1 = fopen ('Euler.dat','w'); % Open output file
for k=1:n % Write to file
  fprintf(f1,fmt_string,x(k),real(z(k)),imag(z(k)));
end
fclose(f1); % Close output file

if graph
  plot(x,real(z),x,imag(z)); % Simple plot in graphical window
  hold on;                   % Wait for closing graphical window
end

printtable([x' real(z)' imag(z)'],'precision',3,'normalized',true)

% COMMENTS

% "Matlab/Octave output precision" should not be confused
% with the number of SIGNIFICANT DIGITS.

% To correctly manage significant digits as required in physics and engineering, use:

% (1) SCIENTIFIC NOTATION 

% fprintf (' %#w.mG ',y);  % where m   is the number of significant digits

% OR

% (2) NORMALIZED SCIENTIFIC NOTATION (always one non-zero digit before decimal point)

% fprintf (' %w.mE ',y);   % where m+1 is the number of significant digits.
% In this case, m is the "Octave/Matlab output precision",
% i.e. the number of digits after the radix (decimal point).

% In both cases (1) and (2): 
% y = real number to be printed respecting significant digit rule.
% w = integer (units of character count) = width of field to print y.
% Typing g instead of G or e instead  E in the format definition toggles
%   the exponent flag respectively in lower or upper case (e or E).

% Both above formats (1) and (2) force to print trailing zeros
% that are significant digits and are useful to perform alignment of columns of numbers.

% MANDATORY : w > m+7 because :
%
% 1 character  for the sign of y
% 1 character  for the decimal point  
% 1 character  for the exponent flag (e or E) 
% 1 character  for the exponent sign 
% 3 characters for the exponent (integer number)

% A good choice is w>m+8 to warrant at least one blank character
% between two numbers on a same line.
