function [n]=columns(A)
    
%% Return number of columns n of a 2 dimension array A(m,n)
%  This shortcut of the size function for 2 dimension array 
%  exists in GNU Octave but not in Matlab.
%  Place this function in search path to enable its use in Matlab.
    
%  See also: rows.m
    
    if(ndims(A)==2)
        n=size(A,2);
    else
        error(['Function columns not applicable if number of ' ...
               'dimensions of argument array ~= 2'])
    end
end % End of function columns