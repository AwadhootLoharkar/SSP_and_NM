function [gapmin,gapmax]=findgaps(w,varargin)

%% FIND GAPS IN QUASIPARTICLE DISPERSION RELATION (BANDS)

%% INPUT

%  w(1:nq,1:L) = array of quasiparticles eigenvalues such that:
%       w(k,i) = eigenvalue of band i at wavevector q(1:3,k)
%            L = number of bands to be plotted

%% OUTPUTS
    
%  gapmin(:) = array of minima of detected gaps (same unit as w)
%  gapmax(:) = array of maxima of detected gaps (same unit as w)
    
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 
    
% if varargin{k}= 'GapTol', 
%                 then varargin{k+1} = real scalar > 0 
%                 = tolerance to recognize a gap (same unit as w)
    
% if varargin{k}= 'RefLevel', 
%                 then varargin{k+1} = real scalar = reflvl
%                 If non-empty, shift the reference level (zero) 
%                 of the ordinates to reflvl
    
% if varargin{k}= 'tol', 
%                 then varargin{k+1} = real scalar > 0 
%                 = tolerance to recognize zero.
%                 In practice, two floating-point numbers
%                 will be considered equal if the absolute value 
%                 of their difference is < tol.
    
% if varargin{k}= 'ylim', 
%                 then varargin{k+1 } = yrange(1:2), real array
%      yrange(1)= lower bound of search interval (same unit as w)
%      yrange(2)= upper bound of search interval (same unit as w)
    
% if varargin{k}= 'Strict_ylim', 
%                 then varargin{k+1 } = Boolean (default = false).
%                 If true, option 'ylim' is applied strictly,
%                 thereby meaning that gaps featuring only gap_min 
%                 or gap_max inside the search interval will not be
%                 part of the output. 

%% INPUT PARSING

pa=inputParser; pa.CaseSensitive = false; % Initialize input parser
addRequired(pa,'w',          @(x) isreal(x) && ismatrix(x));
addParameter(pa,'gaptol',[],                @(x) isreal(x) && x>0);
addParameter(pa,'reflevel',[],  @(x) isscalar(x) && isreal(x));
addParameter(pa,'tol',1e-12,    @(x) isscalar(x) && isreal(x) && x>0);
addParameter(pa,'ylim',[], @(x) isreal(x) && all(size(x)==[1 2]));
addParameter(pa,'Strict_ylim',false, @(x) isscalar(x) && islogical(x));

% Read the argument list knowing the rules
parse(pa,w,varargin{:}); 
w=pa.Results.w; 
gaptol=pa.Results.gaptol; 
reflvl=pa.Results.reflevel; 
tol=pa.Results.tol;
yrange=pa.Results.ylim; 

nq=size(w,1);  % Number of wavevectors
L=size(w,2);   % Number of bands

%% SHIFT REFERENCE LEVEL

if ~isempty(reflvl)
    w=w-reflvl;
end

%% GAPS

if isempty(gaptol)
    amply=abs(min(min(w))-max(max(w)));
    gaptol=10^(floor(log10(amply))-4);
end

fprintf(['\nThreshold to recognize a gap in findgaps.m = %d\n'...
         '(Adjustable through option ''gaptol'').\n'],gaptol)

for i=1:L
    bandmin(i)=min(w(:,i));
    bandmax(i)=max(w(:,i));
end
gapmin=[];gapmax=[];
for i=1:L-1
    if bandmax(i)+gaptol < bandmin(i+1)
        gapmin = [ gapmin bandmax(i)];
        gapmax = [ gapmax bandmin(i+1)];
    end
end
gapmin(abs(gapmin)<tol)=0; % Filtering to set nice zeroes
gapmax(abs(gapmax)<tol)=0;


printtable([ (1:length(gapmin))' gapmin' gapmax'],'Integer',1,...
           'Colname',{'Gap_nb','Gap_min','Gap_max'},...
           'Title','Gap(s) found by findgaps.m:')

%% OPTIONAL: KEEP ONLY GAPS INSIDE [yrange(1),yrange(2)]

if ~isempty(yrange)
    yrange=sort(yrange);
    keptmin=[];keptmax=[];
    for i=1:length(gapmin)
        if gapmin(i) >= yrange(1) && gapmax <= yrange(2)
            keptmin=[keptmin gapmin(i)];
            keptmax=[keptmax gapmax(i)];
        end
        if ~pa.Results.Strict_ylim
            if gapmin(i) <= yrange(1) && gapmax <= yrange(2)
                warning(['Gap with Gap_min outside search interval added to output. ' ...
                         'Toggle option ''Strict_ylim'' to true in order to ' ...
                         'exclude the associated gap from the output.'])
                keptmin=[keptmin gapmin(i)];
                keptmax=[keptmax gapmax(i)];
            end
            if gapmin(i) >= yrange(1) && gapmax >= yrange(2)
                warning(['Gap with Gap_max outside search interval added to output. ' ...
                         'Toggle option ''Strict_ylim'' to true in order to ' ...
                         'exclude the associated gap from the output.'])
                keptmin=[keptmin gapmin(i)];
                keptmax=[keptmax gapmax(i)];
            end
        end
    end
    gapmin=[];      gapmax=[]; % Flush
    gapmin=keptmin; gapmax=keptmax;
end

end % end of function findgaps