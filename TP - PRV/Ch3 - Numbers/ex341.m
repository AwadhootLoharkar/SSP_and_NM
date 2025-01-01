%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%    Num Met 4 Phys - Ex 3.4.1  %
%            1 Oct 24           %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../../start-kit-student")
addpath("../bin")

mystartdefaults;
tol = 1E-12;


%% Question 1
fprintf("\n*** QUESTION 1 ***\n");

PiS = single(pi); % Single precision
PiD = double(pi); % Double precision

fprintf("*** USING %%1.50G ***\n");

fprintf("True value of Pi\n3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282\n");
fprintf("Single precision: 7 (sometimes 8) first digits are correct\n");
fprintf("%1.50E\n", PiS);
fprintf("Double precision: 15 (sometimes 16) first digits are correct\n");
fprintf("%1.50E\n", PiD);

fprintf("\n");
fprintf("*** USING %%1.50E ***\n");

fprintf("True value of Pi\n3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282\n");
fprintf("Single precision: 7 first digits are correct\n");
fprintf("%1.50G\n", PiS);
fprintf("Double precision: 16 first digits are correct\n");
fprintf("%1.50G\n", PiD);


%% Question 2
fprintf("\n*** QUESTION 2 ***\n");

e = single(1); p = 0;
while e ~= 0
    p = p+1;
    e = single(2^(-1*p));
end
eminS = 2^(-1*(p-1));
fprintf("Computed value of eps_0 in simple precision\n")
fprintf("%10.50G\n", eminS);

e = double(1); p = 0;
while e ~= 0
    p = p+1;
    e = double(2^(-1*p));
end
eminD = 2^(-1*(p-1));
fprintf("Computed value of eps_0 in double precision\n")
fprintf("%10.50G\n", eminD);

fprintf("\n");

fprintf("True value of eps_0 in single precision\n")
fprintf("%10.50G\n", eps(single(1.0)));

fprintf("True value of eps_0 in double precision\n")
fprintf("%10.50G\n", eps(double(1.0)));


%% Question 3
fprintf("\n*** QUESTION 3 ***\n");

fprintf("\n\n-----------------------------------------------------------------------------------------------------------\n");
fprintf("| Precision |     n_min     |    n_max     |    eps_0     |       v       |       u       |     e_min     |\n")
fprintf("----------------------------------------------------------------------------------------------------------\n");
nmin = -126; nmax = 127; epsilon = eps(single(1.0)); v = (2-epsilon)*2^nmax; u = 2^nmin; 
fprintf("|  Single   | %e | %e | %e | %e  | %e  | %e  |\n", nmin, nmax, epsilon, v, u, eminS);
fprintf("----------------------------------------------------------------------------------------------------------\n");
nmin = -1022; nmax = 1023; epsilon = eps(double(1.0)); v = (2-epsilon)*2^nmax; u = 2^nmin; 
fprintf("|  Double   | %e | %e | %e | %e | %e | %e |\n", nmin, nmax, epsilon, v, u, eminD);
fprintf("-----------------------------------------------------------------------------------------------------------\n");


%% Question 4
fprintf("\n*** QUESTION 4 ***\n");

fprintf("u-realmin: %e\n", u-realmin());
fprintf("v-realmax: %e\n", v-realmax());
fprintf("\nu = realmin and v = realmax...\n");


%% Question 5
fprintf("\n*** QUESTION 5 ***\n");

fprintf("|  k |      u/2^k      | Sig digits |                             Bit pattern                          |\n");
for k=0:53
    res = u/2^k;
    signDigits=floor(log10(2^(52-k))); % Significant digits
    if (signDigits < 0)
        signDigits = 0;
    end
    bitpattern = float2bin(res);
    fprintf("| %2d | %15e |     %2d     | %s |\n", k, res, signDigits, bitpattern);
end


%% Question 6
fprintf("\n*** QUESTION 6 ***\n");

for k=-308:1:308
    x = 10^k;
    epsilon(1) = 2^(floor(log2(x))-52);
    epsilon(2) = eps(x);
    epsilon(3) = abs(x)*eps; % Upper boundary
end


%% Question 7

%% Question 9

% If the argument of the exponential is more than 709, it will overflow!

%% Question 10

% 1/u create a subnormal number (quite ok)
% Offset binary is inconsistant
