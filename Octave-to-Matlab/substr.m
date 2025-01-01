function [result]=substr(str,m,n)

%% MATLAB TRANSLATION OF OCTAVE FUNCTION substr
%  Extract the substring from str that occurs between the positions m
%  and n, including the characters at those positions.

%%  INPUT
%  str = character string
%  m = integer > 1 and <= length(str)
%  n = integer > m and <= length(str)

%% OUTPUT 
%  result = character string extracted from input str between
%  positions m and n included

result=string(extractBetween(str,m,n)); % Caution: extractBetween
                                        % returns a cell array that
                                        % must be converted into a string

end % End of function substr