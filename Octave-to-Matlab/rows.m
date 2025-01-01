function [m]=rows(A) 

%% Return number of rows m of a 2 dimension array A(m,n)
%  This shortcut of the size function for 2 dimension array 
%  exists in GNU Octave but not in Matlab.
%  Place this function in search path to enable its use in Matlab.
    
%  See also: columns.m

    if(ndims(A)==2)
        m=size(A,1);
    else
        error(['Function rows not applicable if number of ' ...
               'dimensions of argument array ~= 2'])
    end
end % End of function rows