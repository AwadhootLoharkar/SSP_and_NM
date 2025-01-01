function [qpath,bandlabels]=BZpath(step,qs,qe,qs_str,qe_str,ticksep,varargin)

%% GENERATE PATH IN BRILLOUIN ZONE (BZ) 
%% USED AS ABSCISSAS IN A DISPERSION RELATION OF QUASIPARTICLES
%  [-] In 2*pi/spacing units (where spacing = lattice spacing);
%  [-] According to traditional solid state physics representation;
%  [-] Providing a fictitious linearly increasing coordinate
%      for the purpose of a later graphical representation
%      of the dispersion relations of quasiparticles by plotbands.m;
%  [-] Providing abscissas of ticks to be used by plotbands.m.;
%  [-] Providing abscissa tick labels and their positions 
%      to be used by plotbands.m.

%% INPUT
%  step = Step along path in BZ
%  qs(1:3,i) = Start point of segment i 
%  qe(1:3,i) = End   point of segment i 
%  qs_str{i} = Label of start point of segment i 
%  qe_str{i} = Label of end   point of segment i 
%  ticksep(i)= Distance between tick marks along segment i 
%              such that 0 < ticksep(i) <= 1
%  where the index i runs over the number of segments (s).
 
%% OUTPUT
%  qpath(1:3,k) = Cartesian coordinates of wavevector k
%  qpath(4,k)   = Cartesian norm of wavevector k 
%  qpath(5,k)   = Fictitious stricly increasing coordinate 
%                 starting with q(5,1)=0
%                 for the purpose of graphical representation 
%  qpath(6,k)   = Segment identifier of vector qpath(1:3,k)
    
%  bandlabels   = Structure carrying data needed to manage abscissa
%                 labels of a forthcoming plot of dispersion relations
%                 of quasiparticles by plotbands.m
%                 (with s = number of segments of the path in BZ):
    
%  bandlabels.til{1:s}   = string array = tick labels along abscissa axis. 
%  bandlabels.tix(1:s)   = tick positions of abscissa labels
%                          across segments along abscissa axis.
%  bandlabels.qu_str{1:s}= string array = directions of segments according 
%                          to the notation of crystallography (e.g. [110])
%  bandlabels.qu_x(1:s)  = x-positions where character strings
%                          bandlabels.qu_str{1:s} will be centered in plot.
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 
%
% if varargin{k} = 'PlotPath', 
%                  then varargin{k+1} = Boolean  
%                  If true, plot the path segments in 3D

% if varargin{k} = 'File', 
%                  then varargin{k+1} = Character string 
%                  = name of output file
                                   
% if varargin{k} = 'Colors', 
%                  then varargin{k+1 } = 7x3 matrix
%                  = Color order of successive graphical objects
    
% if varargin{k} = 'ShowLegend', 
%                  then varargin{k+1} = Boolean 
%                  If true, show legend  

%% INPUT PARSING

pa=inputParser; pa.CaseSensitive = false; % Initialize input parser

addRequired(pa,'step', @(x) isreal(x) && x>0);
addRequired(pa,'qs',   @(x) isreal(x) && ismatrix(x) && size(x,1)==3);
addRequired(pa,'qe',   @(x) isreal(x) && ismatrix(x) && size(x,1)==3);

addRequired(pa,'qs_str', @(x) iscell(x));
addRequired(pa,'qe_str', @(x) iscell(x));

addParameter(pa,'plotpath',false,  @(x) islogical(x) && isscalar(x));
addParameter(pa,'showlegend',true, @(x) islogical(x) && isscalar(x));
addParameter(pa,'file',[],         @(x) ischar(x));
addParameter(pa,'colors',[],       @(x) isreal(x) && all(size(x)==[7 3]));

parse(pa,step,qs,qe,qs_str,qe_str,varargin{:}); % Read the argument list
                                                % knowing the rules
     step=pa.Results.step;
       qs=pa.Results.qs;                 qe=pa.Results.qe;
   qs_str=pa.Results.qs_str;         qe_str=pa.Results.qe_str;
plot_path=pa.Results.plotpath;  show_legend=pa.Results.showlegend;
     file=pa.Results.file;     color_scheme=pa.Results.colors;

%% CORE JOB

%% Check consistency

[~,s] = size(qs); [~,t] = size(qe); 
if ( s ~= t  )
    error('Numbers of start and end points do not match.\n');
else
    fprintf('\nFunction BZpath: %d segments\n',s);
end

if ( length(ticksep)~=s )
    error(['Length of array ticksep does not match the number of ' ...
           'segments.\n']);
end
if ( any(ticksep<=0) || any(ticksep>1) )
    error('Values of array ticksep should be in ]0,1].\n');
end

%% Initialize tick marks

%  At the origin of abscissas:
bandlabels.tix(1)=0; bandlabels.til{1}=qs_str{1}; l=1; 

%% Initialization of BZ path computation

qu=qe-qs;     % Direction vectors
qu=sign(qu);  % Traditional solid state physics
	      % cristallography vectors along each segment

for is=1:s
    fprintf('Segment %d: direction [%d %d %d]''\n',is,qu(1:3,is))
end

q=zeros(5,1); % Initialize q vector  
p=zeros(5,1); % Initialize vector preceding q
iq=0;         % Initialize number of points for the path in BZ  

%% Loop overs segments

for is=1:s    
    ip=0; cond=0; 
    while cond ~= 1
        xsi=ip*step;
        buf=xsi*qu(:,is);
        % cond=1 if going beyond end point
        % NB: Using norm is not appropriate when path
        %     is such that norm(qe) < norm(qs).
        %     => Use a test that must be true 
        %     for each Cartesian component.
        cond=all(abs(buf)>=abs( (qe(:,is)-qs(:,is)) ) ); 
        
        iq=iq+1;
        if cond==1                   % Last point of segment
            q=qe(1:3,is); q(4)=sqrt(q'*q);
            gap=max(abs(qe(1:3,is)-p(1:3)));
            q(5)=p(5)+gap;
            
            %$$$ Manage labels when last point of segment is found
            ref=bandlabels.tix(l);
            bandlabels.qu_x(is)=ref+0.3*(q(5)-ref);
            bandlabels.qu_str{is}=sprintf('[%d,%d,%d]',...
                                          qu(1,is),qu(2,is),qu(3,is));
            for k=1:round(1/ticksep(is))-1
                l=l+1; bandlabels.tix(l)=ref+k*ticksep(is)*(q(5)-ref);
                if all(qu(:,is)<=0)
                    bandlabels.til{l}=sprintf('%0.2g',1-k*ticksep(is));
                else
                    bandlabels.til{l}=sprintf('%0.2g',k*ticksep(is));
                end
            end
            l=l+1; bandlabels.tix(l)=q(5); bandlabels.til{l}=qe_str{is}; 
            %$$$ End of manage labels
            
        else                        % Standard point of segment
            q=buf(1:3)+qs(1:3,is); q(4)=sqrt(q'*q);
            if iq == 1              % First point of path in BZ
                q(5)=0; p(5)=0;
            else
                if ip == 0
                    q(5)=p(5);      % First point of segment
                else
                    q(5)=p(5)+step; % Standard point of segment
                end
            end
        end
        p=q; ip=ip+1;
        qpath(1:5,iq)=q(1:5); % Store coordinates of new point iq
        qpath(6,iq)=is;       % Segment identifier of point iq
    end % End while
    
end % End of loop over segments

%% Formatted output

if ~isempty(file)
    header{1}='n';
    header{2}='q_1';header{3}='q_2';header{4}='q_3';
    header{5}='q_4';header{6}='q_5';header{7}='q_6';
    printtable([(1:iq)' qpath'],'Integer',[1 7],...
               'ColName',header,'file',file)
end

%% OPTIONAL PLOT

if plot_path
    
    plotopt=plot_preamble; k=plotopt.k; Legend=plotopt.Legend;
    
    % Plot segments of path
    for i=1:s
        k=k+1; Legend{k}=sprintf('segment %d',i);
        pt(1:3,1)=qs(1:3,i); pt(1:3,2)=qe(1:3,i);
        plot3(pt(1,:),pt(2,:),pt(3,:),...
              '-','color',plotopt.color_scheme(i,1:3),...
              'linewidth',3,...
              'DisplayName',Legend{k}); hold on; 
    end
    
    % Plot start points of segments
    for i=1:s
        k=k+1; Legend{k}=''; % Skip this item in Legend
        plot3(qs(1,i),qs(2,i),qs(3,i),...
              'o','color','k',...
              'MarkerSize',plotopt.mks,...
              'MarkerFaceColor','k'); hold on; 
        
        if (strcmp(qs_str{i},'\Gamma'))
            offset=-0.075*max(max(qs));
            text(offset,offset,-offset,qs_str{i},...
                 'fontsize',20,'color','k');
        else
            offset=1.1;
            pos(:)=qs(:,i)*offset; pos(3)=pos(3)+0.2;
            for icar=1:3
                if abs(pos(icar))< 0.05
                    pos(icar)=0.05;
                end
            end
            text(pos(1),pos(2),pos(3),qs_str{i},...
                 'fontsize',20,'color','k');
        end
    end
    
    % Plot end points of segments
    for i=1:s
        k=k+1; Legend{k}=''; % Skip this item in Legend
        plot3(qe(1,i),qe(2,i),qe(3,i),...
              'o','color','k',...
              'MarkerSize',plotopt.mks,...
              'MarkerFaceColor','k'); hold on; 
        
        if (strcmp(qe_str{i},'\Gamma'))
            offset=-0.075*max(max(qs));
            text(offset,offset,-offset,qe_str{i},...
                 'fontsize',20,'color','k');
        else
            offset=1.1;
            pos(:)=qe(:,i)*offset; pos(3)=pos(3)+0.2;
            for icar=1:3
                if abs(pos(icar))< 0.05
                    pos(icar)=0.05;
                end
            end
            text(pos(1),pos(2),pos(3),qe_str{i},...
                 'fontsize',20,'color','k');
        end
    end
        
    % Final settings
    
    axis equal; grid on; legend(Legend,'Location','EastOutside');
    xlabel('\itq\rm_1'); ylabel('\itq\rm_2'); zlabel('\itq\rm_3');
   
    if(~show_legend)
        legend(gca,'off');
    end
    
    view([1,0.5,0.3]); % Point of view
    
end

end  % End of function BZpath
