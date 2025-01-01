%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%   Num Met 4 Phys - Ex 2.12.3  %
%           17 Sept 24          %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../../start-kit-student")

mystartdefaults;

%% Ideal values

func = @(x) exp(-0.2*x).*cos(x);
x = 0:0.1:10; n = length(x);
y = func(x);

%% Adding some noise

noise = rand(1, n)-0.5; % Between -.5 and .5
A = 0.5; % Amplitude of noise
yref = y+A*noise;

%% Plot

figure; hold on;
plot(x, y);
plot(x, yref);

%% Characterize errors

stats = gof(y, yref)