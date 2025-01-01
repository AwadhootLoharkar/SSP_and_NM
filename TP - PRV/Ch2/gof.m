%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%   Num Met 4 Phys - Ex 2.12.3  %
%           18 Sept 24          %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [z] =  gof(y, yref, varargin)
    %% Goodness Of Fit
    % INPUTS
    %   y       array of approximated values
    %   yref    array of reference values
    %
    % OUTPUTS
    %   z       structure array containing the folowing fields
    %     MEAN  Mean of reference values
    %      TSS  Total Sum of Squares
    %      SSE  Sum of Squares of Errors
    %      SSR  Sum of Squares of the Regression
    %       SZ  Sum expected to be Zero if least squares optimization was performed
    %     NSSR  Normalized Sum of Squares of the Regression = SSR/TSS
    %      NSZ  Normalized SZ = SZ/TSS
    %      MSE  Mean of Squares of Errors
    %    NRMSE  Normalized Root Mean of Squares of Errors = sqrt(NMSE)
    %     NMSE  Normalized Mean of Squares of Errors
    %  Rsquare  R-squared or R^2 (coeff. of determination)
    %      FOM  Figure of merit of the optimization
    %


    %% Input parsing
    p = inputParser; p.CaseSensitive = false;
    addRequired(p, 'y',                 @(x) isreal(x) && isvector(x));
    addRequired(p, 'yref',              @(x) isreal(x) && isvector(x));
    addOptional(p, 'tol', 1e-12,        @(x) isreal(x) && isscalar(x) && x>0);
    addOptional(p, 'short', false,      @(x) islogical(x) && isscalar(x));
    parse(p, y, yref, varargin{:});
    y = p.Results.y; yref = p.Results.yref; tol = p.Results.tol; short = p.Results.short;

    if(length(y) - length(yref) ~= 0)
        error("gof: y and yref must be of same length");
    end
    n = length(y);

    %% Computing                        y_i = yref ; \tilde y = y

    MEAN = sum(yref)/n;                 % \bar y
    TSS = sum((yref-MEAN).^2);          % S_T
    SSE = sum((yref-y).^2);             % S_E
    SSR = sum((y-MEAN).^2);             % S_R
    SZ = 2*sum((yref-y).*(y-MEAN));     % S_Z
    NSSR = SSR/TSS;
    NSZ = SZ/TSS;
    MSE = SSE/n;
    NMSE = SSE/TSS;                     % N_E
    NRMSE = sqrt(NMSE);                 % R_E
    RSQUARE = 1 - NMSE;                 % R^2
    FOM = heaviside(RSQUARE)*(1-SZ/SSR);

    %% Checks

    if(SSE+SSR+SZ-TSS > tol)
        fprintf("WARN: S_T is not equal to S_E + S_R + S_Z.\n");
    end

    FOM(abs(FOM)<tol) = 0;

    %% Returned struct

    if(~short) % Full data
        z.MEAN = MEAN;          z.TSS = TSS;        z.SSE = SSE;
        z.SSR = SSR;            z.SZ = SZ;          z.NSSR = NSSR;
        z.NSZ = NSZ;            z.MSE = MSE;        z.NRMSE = NRMSE;
    end

    z.NMSE = NMSE;      z.RSQUARE = RSQUARE;        z.FOM = FOM;
end