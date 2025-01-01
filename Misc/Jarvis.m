function [vertex,varargout]=Jarvis(R,varargin)

%% JARVIS' MARCH:
%% ALGORITHM FINDING THE VERTICES
%% OF THE CONVEX HULL OF A SET OF POINTS IN A PLANE
   
% Ref: R.A. Jarvis, "On the identification of the convex hull 
%      of a finite set of points in the plane". 
%      Information Processing Letters vol 2, 18–21 (1973). 
%      doi:10.1016/0020-0190(73)90020-3
    
%% DEPENDENCE ON EXTERNAL FUNCTIONS
%  JarvisTurn.m
    
%% INPUT
%  
%  R = coordinates of n points in plane x-y  [length]
%      (2Xn array = column order) 
% 
%  R(1,i) = x-coordinate of point i (i=1,...,n)
%  R(2,i) = y-coordinate of point i (i=1,...,n)
    
%% OUTPUT
%
%  vertex = column ordered coordinates of vertices of the 
%           convex hull containing all input points with 
%           the properties enabling plotting a 
%           closed hull as a polygon, i.e:
%             - the vertices are ordered anticlockwise;
%             - the last vertex = first vertex.  
%               Therefore, the correct number of vertices 
%               of the convex hull polygon = columns(vertex)-1
%
%           vertex(1,i) = x-coordinate of vertex i  
%           vertex(2,i) = y-coordinate of vertex i  
    
%  varargout{1} = area of the convex hull polygon [length^2]
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 
    
% if varargin{k} = 'tol', 
%                  then varargin{k+1} = real scalar > 0 
%                  = tolerance to recognize zero.
%                  In practice, two floating-point numbers
%                  will be considered equal if the absolute value 
%                  of their difference is < tol     
       
% if varargin{k} = 'PlotHull',
%                  then varargin{k+1} = Boolean 
%                  If true, plot convex hull
    
% if varargin{k} = 'PlotArrows',
%                  then varargin{k+1} = Boolean 
%                  If true and if PlotHull=true, add arrows showing 
%                  the order of the march on the convex hull polygon
    
% if varargin{k} = 'PlotPoints',
%                  then varargin{k+1} = Boolean  
%                  If true, plot points 
  
% if varargin{k} = 'ColorHull',
%                  then varargin{k+1} = 1x3 real vector
%                  = RGB color of hull and arrows
    
% if varargin{k} = 'ColorPoints',
%                  then varargin{k+1} = 1x3 real vector
%                  = RGB color of point markers
    
% if varargin{k} = 'ColorStartPoint',
%                  then varargin{k+1} = 1x3 real vector
%                  = RGB color of first vertex of Jarvis' march
    
% if varargin{k} = 'Verbose',
%                  then varargin{k+1} = Boolean 
%                  If true, display of some computation details
    
%% INPUT PARSING

p=inputParser; p.CaseSensitive = false; % Initialize input parser

addRequired(p,'R', ...
             @(x) isreal(x) && ismatrix(x) && size(x,1)==2);

addParameter(p,'tol',1e-12,...
             @(x) isreal(x) && isscalar(x) && x>0);
addParameter(p,'verbose',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plotpoints',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plothull',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plotarrows',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'colorpoints',[0 0 1],...
             @(x) isreal(x) && all(size(x)==[1 3])...
                  && all(x>=0) && all(x<=1) );
addParameter(p,'colorhull',[1 0 0],...
             @(x) isreal(x) && all(size(x)==[1 3])...
                  && all(x>=0) && all(x<=1) );
addParameter(p,'colorstartpoint',[0 1 0],...
             @(x) isreal(x) && all(size(x)==[1 3])...
                  && all(x>=0) && all(x<=1) );

parse(p,R,varargin{:}); % Read the argument list
R=p.Results.R;
tol=p.Results.tol;
verbose=p.Results.verbose;
plot_points=p.Results.plotpoints;
plot_hull=p.Results.plothull;
plot_arrows=p.Results.plotarrows;
color_scheme(1,1:3)=p.Results.colorpoints;
color_scheme(2,1:3)=p.Results.colorhull;
color_scheme(3,1:3)=p.Results.colorstartpoint;

%% PARSING varargout

if max(nargout,1) - 1 > 1;
    error(['Function Jarvis: can not deal '...
           'with more than 1 optional output field.'])
end

%% CORE JOB
 
[~,n]=size(R);

if (n<3)
    error(['Function Jarvis: there must be a least 3 ' ...
           'points.\n']);
end
if (plot_points || plot_hull)
    S=R; % Need a backup of R for plotting
end

[ymin,lowest] = min(R(2,:)); % Ordinate and index of lowest point

vertex=R(:,lowest); % Lowest point = first vertex

R = [R(:,1:lowest-1),R(:,lowest+1:n)]; % Remove first vertex from
                                       % list where to search the
                                       % second vertex

step=1;  a=vertex;  target=vertex; % Initialization
unclosed_loop=true;

while (unclosed_loop); % Looping until back to starting point
    
    b=R(:,1); k=1; % Start with first point of the updated search list
    
    [~,n]=size(R);
    for i=2:n % Search the point that is most on the right
              % relatively to current direction
        c=R(:,i);
        turn=JarvisTurn(a,b,c);
        switch turn
        % case -2 % Points b and c colinear in opposite
                  % directions can not happen in this algorithm 
                  % => No action.
          case -1 % If point c is on the right of point b,
            buf=b; b=c; c=buf; % permute b and c
            k=i; 
          case 0  % Identical points are duplicated in list of vertices
            buf=b; b=c; c=buf; % permute b and c
            k=i; 
        % case 1  % Keep b as possible vertex until end of for <=> No action.
          case 2  % Points b and c colinear in same direction
            if (norm(b-a) <= norm(c-a)) % If c most distant from a,
                buf=b; b=c; c=buf;      % permute b and c.
                k=i;
            end  
        end     
    end
    vertex=[vertex b]; % Next vertex found = last value of b
                       % = point most on the right 
                       %   relatively to current direction
   
    a=b; % Last found vertex becomes the next pivot 
    
    if(norm(a-target) < tol)
        unclosed_loop=false; % Closed loop => end of while
    end
    
    R(:,k)=[]; % Remove the last found vertex from the search list
    if(step == 1)            % After finding the second vertex,
        R=[ vertex(:,1) R ]; % reintroduce first vertex as first
                             % element of the search
                             % list to enable to find it again 
                             % as target terminating the while loop.
    end
    step=step+1; % Increment the number of hull vertices
    
end % End of while (unclosed_loop)

%% AREA OF CONVEX HULL POLYGON 
%  by anticlockwise decomposition in signed triangular areas

%  Basic principle:
%  A_i   = det([vertex(:,i) vertex(:,i+1)]) = cross-product
%        = signed area of parallelogram spanned
%          by vectors vertex(:,i) and vertex(:,i+1)
%        = signed area inside parallelogram whose vertices are:
%          origin, vertex(:,i),vertex(:,i+1)
%          and vertex(:,i) + vertex(:,i+1)
%  A_i/2 = signed area of triangle defined by the 3 points
%          origin, vertex(:,i), vertex(:,i+1),

%  NB: This method is implemented in Matlab/Octave as the 
%      built-in function polyarea invoked as follows:
%      polyarea(vertex(1,:),vertex(2,:));

A=det([vertex(:,step) vertex(:,1)]); % Last triangle
for i=1:step-1
    A = A + det([vertex(:,i) vertex(:,i+1)]);
end
varargout{1}=A/2;

if (verbose)
    % The number of vertices = step-1 
    % because the last vertex = first vertex
    fprintf('Jarvis: Convex hull = %d vertex polygon.\n',step-1)
    fprintf('        Area inside convex hull = %d\n',varargout{1})
end

%% OPTIONAL PLOTS

if (plot_points || plot_hull)
     
    h = findall(0,'type','figure');
    if isempty(h)    
        figure('NumberTitle', 'off',...
               'name',['Anticlockwise Jarvis march  '...
                       'starting at filled circle'],...
               'defaultAxesFontName','Cambria',...
               'defaultAxesFontSize',18)
    end
    
    if ismatlabshell
        lw=2.5; mhs=1.1; mks=6;   mid=0.625;
    end
    if isoctaveshell
        lw=1.5; mhs=0.1; mks=3.5; mid=0.575;
    end
    
    if (plot_hull)
        plot(vertex(1,:),vertex(2,:),'-',...
             'Color',color_scheme(1,1:3),...
             'Linewidth',lw); hold on;
        if (plot_arrows)
            for i=1:step-1
                quiver(vertex(1,i),vertex(2,i),...
                       (vertex(1,i+1)-vertex(1,i))*mid,...
                       (vertex(2,i+1)-vertex(2,i))*mid,...
                       'Color',color_scheme(1,1:3),...
                       'MaxHeadSize',mhs,...
                       'Linewidth',lw,...
                       'AutoScaleFactor',1);
            end
        end
    end
    
    if(plot_points)
        % Use backup S for plotting because
        % R was mofified by algorithm
        plot(S(1,:),S(2,:),'o',...
             'Color',color_scheme(2,1:3),...
             'MarkerSize',mks); hold on;
        % First vertex somewhat larger and in another color
        plot(vertex(1,1),vertex(2,1),'o',...
             'Color',color_scheme(3,1:3),...
             'MarkerSize',1.5*mks,...
             'MarkerFaceColor',color_scheme(3,1:3)); hold on;
    end
    
    axis equal; 
    xlabel('\itx\rm_1'); ylabel('\itx\rm_2');
    xlim([floor(min(S(1,:))),ceil(max(S(1,:)))])
    ylim([floor(min(S(2,:))),ceil(max(S(2,:)))])
    
end

end % End of function Jarvis



