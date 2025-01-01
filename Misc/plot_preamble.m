function [plotopt]=plot_preamble()
%% PREAMBLE OF PLOT OPTION(s) CALLED INSIDE A FUNCTION
    
%% Purpose: consistent management of legend items
%% and plot parameters when adding plot objects 
%% by successive calls to various functions
    
%  plotopt is a structure array carrying the following fields:
%      plotopt.k   = current number of items in the legend
%                    of the current plot window;
%      plotopt.Legend = array of legend strings in the current plot window;
%      plotopt.color_scheme = current color scheme (7x3 real array);
%      plotopt.lw  = line width;
%      plotopt.elw = edge line width, typically used in function patch;
%      plotopt.mhs = maxheadsize, used by functions quiver and quiver3;
%      plotopt.mks = markersize.

%  USE: Put this file in your Matlab/Octave path.  Call this preamble
%       in some current_script_or_function.m before the plot commands
%       that are intended to be added (maybe optionally) to some
%       already existing plot window that might have been created
%       outside current_script_or_function.m. The legend items and
%       other plot parameters will then be consistently managed using
%       the data carried by the cell array plotopt.
    
    
%  After calling this function:

%  [1] A figure is opened if none was already open;
%  [2] A legend is initialized as the array plotopt.Legend
%      if none was already activated;
%  [3] The legend of the eventually already open window has been
%      retrieved inside the cell array plotopt.Legend{1:k};
%  [4] plotopt.k is the current number of items in the legend
%      of the current plot window;
%  [5] plotopt.k>=1 because at least a legend title has been created;
%  [6] Color scheme has been retrieved by the command:
%      plotopt.color_scheme=get(gca,'colororder');
%      If plotopt.color_scheme was not yet defined,
%      a color_scheme is defined (starting by blue, red, green);
%  [7] Some default plot parameters have been set differently 
%      according to the current shell (Matlab or Octave), namely:
%      plotopt.lw, plotopt.elw, plotopt.mhs, plotopt.mks
 
%  NB: This preamble was found to ensure the Matlab/Octave
%      compatibility of the management of legend items for both
%      software versions issued before April 15th, 2024.
    
    
h = findall(0,'type','figure');
if isempty(h) 
    % Open a figure if none already exists
    fig=figure('defaultAxesFontSize',18);
    initialize_legend3D; % Initialization of a legend 
end

% A figure is now surely already open, waiting for input
try
    hleg=legend(); % Was a legend defined in the already open figure?              
catch
    initialize_legend3D; % Initialization of a legend 
                         % if not existing in already open figure
    hleg=legend();  
end

if isempty(hleg)
    error('Empty legend not expected.'); % Looks like overkill test.  
                                         % True purpose: checking
                                         % consistency of built-in legend()
                                         % when software version changes.
else
    % Retrieve legend of current figure
    % Matlab storing Legend items as a single row, 
    % Octave storing them as a single column 
    % => max(size(Legend)) works for both Matlab & Octave
    plotopt.Legend = get(legend(),'string'); 
    plotopt.k=max(size(plotopt.Legend)); % Inherited number of items in Legend
end

plotopt.color_scheme=get(gca,'colororder');
if isempty(plotopt.color_scheme)
    plotopt.color_scheme=[ % Start with blue, red, green
                   0.0000   0.4470   0.7410 % (1) Tropical Blue
                   0.6350   0.0780   0.1840 % (2) Sienna
                   0.4660   0.6740   0.1880 % (3) Grass green
                   0.3010   0.7450   0.9330 % (4) Azure blue
                   0.8500   0.3250   0.0980 % (5) Deep orange
                   0.9290   0.6940   0.1250 % (6) Deep yellow
                   0.4940   0.1840   0.5560 % (7) Violet
                   ];
end

% lw = line width, elw = edge line width, 
% mhs = max head size, mks = marker size
if ismatlabshell
    plotopt.lw=3; plotopt.elw=2;   plotopt.mhs=0.5; plotopt.mks=9;
end
if isoctaveshell
    plotopt.lw=2; plotopt.elw=1.5; plotopt.mhs=0.1; plotopt.mks=5;
end

end % end of function plot_preamble




