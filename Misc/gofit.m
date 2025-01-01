
function [z]=gofit(y,yref,varargin)
    
%% GOODNESS OF FIT 
%  Evaluate quality of an array of approximation values 
%  in fitting an array of reference values 
    
%% INPUTS
    
%  y    = array of approximation values
%  yref = array of reference values
    
%% Recognized options in varargin 
% (uppercases for readability are optional)
    
% if varargin{k}= 'tol', 
%                 then varargin{k+1} = real scalar > 0 
%                 = tolerance to recognize zero.
%                 In practice, two floating-point numbers
%                 will be considered equal if the absolute value 
%                 of their difference is < tol
    
% if varargin{k}= 'short', 
%                 then varargin{k+1} = logical (Boolean).
%                 If true, make a shorter output structure array
%                 (Didactical to illustrate the use of varargin).
    
%% OUTPUTS
    
%  z = structure array containing the fields:
%
%          MEAN: Mean of reference values
%           TSS: Total Sum of Squares
%           SSE: Sum of Squares of Errors
%           SSR: Sum of Squares of the Regression
%            SZ: Sum expected to be Zero if least squares
%                optimization was performed 
%          NSSR: Normalized Sum of Squares of the Regression = SSR/TSS
%           NSZ: Normalized SZ = SZ/TSS
%           MSE: Mean of Squares of Errors
%         NRMSE: Normalized Root Mean of Squares of Errors = sqrt(NMSE)
%          NMSE: Normalized Mean of Squares of Errors
%       Rsquare: R-squared or R^2 (coeff. of determination)
%           FOM: Figure of merit of the optimization
    
%% INPUT PARSING

p=inputParser; p.CaseSensitive = false; % Initialize input parser
addRequired(p,'y', ...
            @(x) isreal(x) && isvector(x));
addRequired(p,'yref', ...
            @(x) isreal(x) && isvector(x));
addParameter(p,'tol',1e-12,...
             @(x) isreal(x) && isscalar(x) && x>0);
addParameter(p,'short',false,...
             @(x) islogical(x) && isscalar(x)); 

parse(p,y,yref,varargin{:}); % Read the argument list
y=p.Results.y;  
yref=p.Results.yref; 
tol=p.Results.tol;
short=p.Results.short;

    
%% WORK: Compute MEAN, TSS, SSE, SSR, SZ, NSSR, NSZ, MSE, NRMSE, NMSE,
n = length(y)
MEAN = mean(yref);
TSS  = sum((yref-MEAN).^2);
SSE  = sum((y-yref).^2);
SSR  = sum((y-yref-MEAN).^2);
SZ   = sum((y-yref).*(y-MEAN));
NSSR = SSR/TSS;
NSZ  = SZ/TSS;
MSE  = SSE/n;
NRMSE= sqrt(MSE);
NMSE = NRMSE/MEAN;
Rsquare = 1 - SSE/TSS;
FOM = heaviside(Rsquare) * (1 - abs(SZ/SSR));



% Fields of the returned structure array

if ~short
    z.MEAN     = MEAN;
    z.TSS      = TSS;
    z.SSE      = SSE;
    z.SSR      = SSR;
    z.SZ       = SZ;
    z.NSSR     = NSSR;
    z.NSZ      = NSZ;
    z.MSE      = SSE/n;
    z.NRMSE    = NRMSE;
end
z.NMSE     = NMSE;
z.Rsquare  = Rsquare;
z.FOM      = FOM;
    
end % End of function gof
