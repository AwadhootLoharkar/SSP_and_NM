function [MATLAB]=ismatlabshell()
    MATLAB=~isempty(strfind(evalc('ver'),'MATLAB')); % True if MATLAB shell
end