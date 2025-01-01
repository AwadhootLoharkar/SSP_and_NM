%% mystartdefaults.m

%  To be invoked as first line of an *.m main script (not inside functions)

scriptinit;        % Clear all traces of previous commands
startup;           % Update specific paths 
                   % (overcaution if the environment did not change)
plotpreferences;   % (Re)-load preferences for plots
physicalconstants; % (Re)-load physical constants
