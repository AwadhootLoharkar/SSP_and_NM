function [boolean]=isintegervalue(n)
    
%%  DETERMINE IF A FLOATING POINT NUMBER MATCHES AN INTEGER VALUE
    
    p=inputParser; p.CaseSensitive = false; % Initialize parser
    addRequired(p,'n', @(x) isreal(x)); % Only real numbers as input
    parse(p,n);   % Read the argument list
    n=p.Results.n;
    
    boolean=~mod(n,1); % True if n is an integer value
    
end % End of function isintegervalue