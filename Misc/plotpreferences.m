%% CONFIGURATION OF PREFERENCES FOR USER INTERFACE AND PLOTS

%  Function "get" and "set" are not case sensitive. 
%  Upper cases for readability.

%  groot is the "handle index" (identifier) of the object 
%  that is "parent" of all plot objects of the session.

%  get(groot,'factory') % List factory-defined plot configuration
%  get(groot,'factoryObjectType') % List factory-defined properties
                                  % of 'ObjectType'
%  Examples of 'ObjectType' are: 
%   'Axes', 'Figure', 'Image', 'Line','Surface', 'Text',
%   'uicontrol' (user interface), etc.

% set(groot,'property',value) has a global effect,
% therefore, only the variables below must be declared as global
global default_colors brg_colors crude_colors

%% USER INTERFACE (UI) DEFAULTS
set(groot,'defaultUicontrolFontSize',13);

%% DEFAULT POSITION AND SIZE OF FIGURE WINDOW
golden_ratio=(1+sqrt(5))/2;
default_screen_size=...
    get(groot,'ScreenSize'); % Get screen dimensions [pixels]
default_fig_pos(4)=...
    default_screen_size(4)/3; % Height [pixels]
default_fig_pos(3)=...
    default_fig_pos(4)*golden_ratio; % Width [pixels]
set(groot,'DefaultFigurePosition',...
          [default_screen_size(3)/2-default_fig_pos(3)/2,...
           default_screen_size(4)/2-default_fig_pos(4)/2,...
           default_fig_pos(3), default_fig_pos(4)]);

%% LOCATION OF ABSCISSA TICKS & LABELS 
%  Possible values:'top','bottom',left','right','origin'.
%  NB:'origin' puts labels of abscissae
%      along a y=0 line inside the plot frame.
set(groot,'defaultAxesXaxisLocation','bottom'); 

%% ADJUST LINE & CHARACTER PROPERTIES 
%% FOR VIDEOPROJECTION IN CLASSROOM

if ismatlabshell
    set(groot,'DefaultAxesFontName','Cambria')
    set(groot,'defaultAxesFontWeight','normal');
    set(groot,'defaultAxesFontSize',18); 
    set(groot,'defaultLegendFontName','Cambria'); 
    set(groot,'defaultLegendFontSize',18); 
    set(groot,'defaultLineLineWidth',2);
    set(groot,'defaultAxesLineWidth',2); 
end

if isoctaveshell
    set(groot,'defaultAxesFontName','Cambria');
    set(groot,'defaultAxesFontWeight','normal');
    set(groot,'defaultAxesFontSize',16); 
    set(groot,'defaultLegendFontName','Cambria'); 
    set(groot,'defaultLegendFontSize',16); 
    set(groot,'defaultLineLineWidth',1.5);
    set(groot,'defaultAxesLineWidth',1.5); 
end

%% FACTORY-DEFINED COLOR ORDER FOR COLORING CURVES
%  Color order of successive curves must be defined 
%  by a 7x3 matrix of values between 0 and 1

%  NB: gca="get current axes" identifier  

%  default_colors = get(gca,'colororder'); % Factory setup
%  NB: Uncommented previous line returns the same definition
%      of the array default_colors as the next command

default_colors=[ % Factory setup
                 0.0000   0.4470   0.7410 % (1) Tropical Blue
                 0.8500   0.3250   0.0980 % (2) Deep orange
                 0.9290   0.6940   0.1250 % (3) Deep yellow
                 0.4940   0.1840   0.5560 % (4) Violet
                 0.4660   0.6740   0.1880 % (5) Grass green
                 0.3010   0.7450   0.9330 % (6) Azure blue
                 0.6350   0.0780   0.1840 % (7) Sienna
               ];

%% USER-DEFINED COLOR ORDERS

brg_colors=[ % Start with blue, red, green
             0.0000   0.4470   0.7410 % (1) Tropical Blue
             0.6350   0.0780   0.1840 % (2) Sienna
             0.4660   0.6740   0.1880 % (3) Grass green
             0.3010   0.7450   0.9330 % (4) Azure blue
             0.8500   0.3250   0.0980 % (5) Deep orange
             0.9290   0.6940   0.1250 % (6) Deep yellow
             0.4940   0.1840   0.5560 % (7) Violet
           ];

crude_colors=[ % Unsoftened colors
               0    0    0      % (1) Black		    
               1    0    0      % (2) Red
               0    0    1      % (3) Blue
               0    0.5  0      % (4) Dark Green
               0.9  0.5  0.1    % (5) Orange
               0    0.75 0.75   % (6) Turquoise
               0.5  0.5  0.5    % (7) Grey
             ];

%% SET DEFAULT COLOR ORDER
set(groot,'defaultAxesColorOrder',brg_colors);

%% ANONYMOUS FUNCTIONS FOR TRACING HORIZONTAL 
%% AND VERTICAL LINES ACROSS THE PLOT FRAME
%  NB: Anonymous functions have a global effect

% yline(yval) traces horizontal line at y(x)=yval using current xlim
yline = @(yval, varargin) line(xlim, [yval yval], varargin{:}); 

% xline(xval) traces vertical   line at    x=xval using current ylim
xline = @(xval, varargin) line([xval xval], ylim, varargin{:}); 