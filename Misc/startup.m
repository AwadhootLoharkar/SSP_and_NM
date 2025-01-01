%% startup.m: MATLAB/OCTAVE STARTUP FILE

% For Matlab: to be stored in MATLAB directory
%             (typically $HOME/Documents/MATLAB)

% For Octave: to be stored as the hidden resource file
%             $HOME/.octaverc

%% User's specific paths

addpath('C:\Users\awadh\github\QSSP & NM');    % Add a directory to search path 

if ismatlabshell
    addpath('C:\Users\awadh\github\QSSP & NM\Octave-to-Matlab') 
    % Path to MATLAB translations of some 
    % GNU OCTAVE functions not built in MATLAB
end

if isoctaveshell
    addpath('C:\Users\awadh\github\QSSP & NM\Matlab-to-Octave')
    % Path to GNU Octave translations of some  
    % MATLAB functions not built in GNU OCTAVE
end