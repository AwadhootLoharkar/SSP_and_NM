%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%     Num Met 4 Phys - Ex 3.8   %
%           15 Oct 24           %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../../start-kit-student")
addpath("../bin")

mystartdefaults;
tol = 1E-12;

%% Question 1+2

A = [1 2 3 4; 8 5 7 2; 9 5 4 6; 0 7 2 6];
[Q, R] = GramSchmidtClassic(A)
[Q, R] = GramSchmidtStable(A)


A-Q*R