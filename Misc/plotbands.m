function plotbands(q,bandlabels,w,varargin)

%% PLOT OF QUASIPARTICLE DISPERSION RELATIONS (BANDS)
%% ALONG SOME PATH IN BRILLOUIN ZONE (BZ)

%% Dependence on external functions
%  [-] BZpath.m defining path in BZ (q) and related labels for
%      plotting (bandlabels);
%  [-] Any script or function computing quasiparticle eigenvalues (w)
%      along the path defined by BZpath..

%% INPUTS DEFINING ABSCISSAS & THEIR LABELS = OUTPUT OF BZpath.m

%% q:
%  The wavevector array returned by BZpath.m 
%  is organised as follows:
%  q(1:3,k) = Cartesian coordinates of wavevector k
%  q(4,k) = Cartesian norm of wavevector k 
%  q(5,k) = Fictitious stricly increasing coordinate 
%           starting with q(5,1)=0
%           for the purpose of graphical representation 
%  q(6,k) = Segment identifier of vector qpath(1:3,k)

%  NB: 
%  [1] The wavevector array q must of course be identical to the one 
%      that has been used to compute the band structure stored in w.
%  [2] q(1:4,:) are not relevant in this module. 
%      They are kept to offer a straightforward use
%      of the array q(1:6,:) output by BZpath.m.

%% bandlabels:
%  Structure returned by BZpath.m 
%  carrying data needed to manage abscissa labels 
%  of the intended plot of dispersion relations of quasiparticles
%  (with s = number of segments of the path in BZ):
    
%  bandlabels.til{1:s}   = string array = tick labels along abscissa axis. 
%  bandlabels.tix(1:s)   = tick positions of abscissa labels
%                          across segments along abscissa axis.
%  bandlabels.qu_str{1:s}= string array = directions of segments according 
%                          to the notation of crystallography (e.g. [110])
%  bandlabels.qu_x(1:s)  = x-positions where character strings
%                          bandlabels.qu_str{1:s} will be centered in plot.


%% INPUT DEFINING THE MULTIVALUED ORDINATES = EIGENVALUES

%  w(1:nq,1:L) = array of quasiparticles eigenvalues such that:
%       w(k,i) = eigenvalue of band i at wavevector q(1:3,k)
%                along the path in BZ defined by BZpath.m
%            L = number of bands to be plotted.
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 

% if varargin{k} = 'Color', 
%                  then varargin{k+1 } = [R G B]
%                  = RGB encoding of the color of curves.
%                  Default=[0 0 0] (black).
    
% if varargin{k}= 'FermiLevel', 
%                 then varargin{k+1} = Fermi_level (real scalar)
%                 If non-empty, add the Fermi level 
%                 as a red dashed line in plot (if in frame).

% if varargin{k} = 'Font', 
%                  then varargin{k+1} = Character string 
%                  = font used for labels. Default = 'Cambria'.
    
% if varargin{k} = 'FontSize', 
%                  then varargin{k+1} = integer
%                  = fontsize used for labels. Default = 12.
%                  Best looking results: between 12 and 15.
    
% if varargin{k}= 'GapTol', 
%                 then varargin{k+1} = real scalar > 0 
%                 = tolerance to recognize a gap (same unit as w).
    
% if varargin{k}= 'RefLevel', 
%                 then varargin{k+1} = reflvl (real scalar)
%                 If non-empty, shift the reference level (zero) 
%                 of the ordinates to reflvl.
    
% if varargin{k}= 'ShowGap', 
%                 then varargin{k+1} = Boolean
%                 If true, draw dotted lines to show gaps.
    
% if varargin{k}= 'tol', 
%                 then varargin{k+1} = tol (real scalar > 0) 
%                 = tolerance to recognize zero.
%                 In practice, two floating-point numbers
%                 will be considered equal if the absolute value 
%                 of their difference is < tol
    
% if varargin{k}= 'VacuumLevel', 
%                 then varargin{k+1} = E_vac (real scalar)
%                 If non-empty, add the vacuum level E_vac 
%                 as a grass-green dashed line in plot (if in frame).
    
% if varargin{k} = 'xunit', 
%                  then varargin{k+1} = Character string 
%                  = label of abscissa units.
%                  Default = '\bfq\rm [2\pi/a]'

% if varargin{k} = 'yunit', 
%                  then varargin{k+1} = Character string 
%                  = label of ordinate units.
%                  Default = '\itE\rm [eV]' 
%                  Typical choice for electrons: '\itE\rm [eV]'
%                  Typical choice for phonons:   '\nu [THz]'
%                                             or '\omega [THz]'
    
% if varargin{k} = 'ylim', 
%                  then varargin{k+1 } = yrange(1:2), real array
%      yrange(1) = min. ordinate in plot (same unit as w)
%      yrange(2) = max. ordinate in plot (same unit as w)

%% INPUT PARSING

pa=inputParser; pa.CaseSensitive = false; % Initialize input parser
addRequired(pa,'q',          @(x) isreal(x) && ismatrix(x));
addRequired(pa,'bandlabels', @(x) isstruct(x));
addRequired(pa,'w',          @(x) isreal(x) && ismatrix(x));

addParameter(pa,'color',[0 0 0], ...
             @(x) isreal(x) && all(size(x)==[1 3]));
addParameter(pa,'font','Cambria',@(x) ischar(x));
addParameter(pa,'fontsize',12',  @(x) isintegervalue(x));
addParameter(pa,'fermilevel',[], @(x) isscalar(x) && isreal(x));
addParameter(pa,'gaptol',[],     @(x) isreal(x) && x>0);
addParameter(pa,'reflevel',[],   @(x) isscalar(x) && isreal(x));
addParameter(pa,'showgap',false, @(x) isscalar(x) && islogical(x));
addParameter(pa,'tol',1e-12,     @(x) isscalar(x) && isreal(x) && x>0);
addParameter(pa,'vacuumlevel',[],@(x) isscalar(x) && isreal(x));
addParameter(pa,'xunit','\bfq\rm [2\pi/a]', @(x) ischar(x));
addParameter(pa,'yunit','\itE\rm [eV]',     @(x) ischar(x));
addParameter(pa,'ylim',[], @(x) isreal(x) && all(size(x)==[1 2]));

% Read the argument list knowing the rules
parse(pa,q,bandlabels,w,varargin{:}); 
q=pa.Results.q;             bandlabels=pa.Results.bandlabels;
w=pa.Results.w;             color=pa.Results.color;
font=pa.Results.font;       fontsize=pa.Results.fontsize;
gaptol=pa.Results.gaptol;   reflvl=pa.Results.reflevel;
showgap=pa.Results.showgap; tol=pa.Results.tol;
xunit=pa.Results.xunit;     yunit=pa.Results.yunit;
yrange=pa.Results.ylim;
E_F=pa.Results.fermilevel;  E_vac=pa.Results.vacuumlevel;

%% CHECK INPUT CONSISTENCY

nq=size(q,2); L=size(w,2);
if nq ~= size(w,1)
    error('BZ path and band structure are incompatible.')
end

%% BANNER

fprintf(['\nFont size effects in plotbands.m:\n'...
         'Resizing the plot window may improve '...
         'how labels are displayed.\n'])

%% SHIFT REFERENCE LEVEL

if ~isempty(reflvl)
    w=w-reflvl;
end

%% PLOT BANDS

bands=figure('NumberTitle', 'off',...
             'defaultAxesFontName',font,...
             'defaultAxesFontSize',fontsize,...
             'defaultAxesFontWeight','normal',...
             'name','Quasiparticle dispersion relations');

plot(q(5,1:nq),w(1:nq,1:L),'color',color); hold on; grid('off');

xlim([ q(5,1), q(5,nq)]); xl = xlim; 
amplx = abs(xl(2)-xl(1)); % Amplitude of abscissas in plot units.
                          % Since abscissas of plot  
                          % always start at 0, 1/amplx is the
                          % scaling factor to normalized units.

if isempty(yrange)
    ylim([ floor(min(min(w))), ceil(max(max(w))) ]);
else
    ylim(sort(yrange));
end

%% PLOT BOUNDARIES OF SEGMENTS
% sbf = Segment Boundary Function
% Trick to allow drawing boundaries of segments of the BZ path 

sbf(:)=1.1*ceil(max(max(abs(w)))) * (-1).^round(q(6,:)); 
plot(q(5,:),sbf,'k','linewidth',3); hold on; % Segment boundaries
                                             % in thicker black lines

%% GAPS
%  This section makes the same job as findgaps.m
%  The call to findgaps.m is avoided so that 
%  plotting does not depend on another external function 
%  processing the same data.

if showgap
    if isempty(gaptol)
         amply=abs(min(min(w))-max(max(w)));
         gaptol=10^(floor(log10(amply))-4);
    end
    
    fprintf(['\nThreshold to recognize a gap in plotbands.m = %d\n'...
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
               'Title','Gap(s) found by plotbands.m:')
    for i=1:length(gapmin)
        line(xlim,[gapmin(i) gapmin(i)],'linestyle',':','color','k'); 
        line(xlim,[gapmax(i) gapmax(i)],'linestyle',':','color','k'); 
    end
end

%% FERMI LEVEL
if ~isempty(E_F)
    line(xlim,[E_F,E_F],'linestyle','--',...
         'color',[0.6350   0.0780   0.1840]); % Sienna
    % Text put to the right of the frame 
    % to avoid overwriting by ticklabels
    text('Position',[xl(2)+amplx/100,E_F],...
         'String','\itE_F',...
         'Color',[0.6350   0.0780   0.1840],...
         'HorizontalAlignment','left',...
         'VerticalAlignment','middle',...
         'FontName',font,...
         'FontSize',fontsize,...
         'FontWeight','normal')
end

%% VACUUM LEVEL
if ~isempty(E_vac)
    line(xlim,[E_vac,E_vac],'linestyle','--',...
         'color',[0.4660   0.6740   0.1880]); % Grass green
    % Text put to the right of the frame 
    % to avoid overwriting by ticklabels
    text('Position',[xl(2)+amplx/100,E_vac],...
         'String','\itE_{vac}',...
         'Color',[0.4660   0.6740   0.1880],...
         'HorizontalAlignment','left',...
         'VerticalAlignment','middle',...
         'FontName',font,...
         'FontSize',fontsize,...
         'FontWeight','normal')
end

%% RESIZE THE FRAME
%  More space is needed below the plot frame to accomodate  
%  an additional line to write the labels 
%  of the directions of each segment.

ax=gca; % get current axes
framepos=get(ax,'Position');  % Position of frame inside window
                              % (window size units \in [0,1])
                              % where framepos(1)=left 
                              %       framepos(2)=bottom 
                              %       framepos(3)=width 
                              %       framepos(4)=height

fraction=0.1; % Parameter suitable for fontsize = 12pt
delta=framepos(4)*fraction;
framepos(2)=framepos(2)+delta;% Increase bottom
framepos(4)=framepos(4)-delta;% Reduce height
set(ax,'Position',framepos);  % Resized frame                          

%% STANDARD ABSCISSA AND ORDINATE LABELS
%  Use normalized units for a rigid vertical position
%  independent of the data units

ax.XLabel.Units='normalized';
ylabel(yunit);
xlh=xlabel(xunit);        % Handle to modify xlabel position below
xlh.Position(2)=-fraction;% Reference = top of string of label
                          % relatively to bottom of frame

%% TICKS ACROSS SEGMENTS AND THEIR LABELS

set(ax,'xtick',bandlabels.tix,'xticklabel',bandlabels.til); 

%% ADDITIONAL LINE OF ABSCISSA LABELS: DIRECTIONS OF SEGMENTS
%  Use normalized units for a rigid vertical position
%  of the additional line of abscissa labels
%  (independent of the data units).
%  Reference = top of string of text relatively to bottom of frame

for i=1:length(bandlabels.qu_x)
    text('Units','Normalized',...
         'Position',[bandlabels.qu_x(i)/amplx,-0.6*fraction],...
         'String',bandlabels.qu_str{i},...
         'HorizontalAlignment','center',...
         'VerticalAlignment','top',...
         'FontName',font,...
         'FontSize',fontsize,...
         'FontWeight','normal')
end
end % end of function plotbands