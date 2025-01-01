function [vertex,edge,face,vn]=WignerSeitz(a,varargin)
    
%% COMPUTE AND PLOT WIGNER-SEITZ PRIMITIVE CELL
%  Wigner-Seitz Cell = primitive cell containing a single lattice node.
%  In reciprocal space, Wigner-Seitz cell = Brillouin Zone.
    
%% Dependence on external functions
%  Jarvis.m
%  JarvisTurn.m
%  lattice.m
%  RotationAxe.m
%  plot_preamble.m
%  initialize_legend3D;
%  For MATLAB: translations of Octave's rows.m and columns.m
%  For OCTAVE: initialize_legend3D.m
    
%% Inputs:
%  a = 3X3 array of lattice basis vectors in column order [length]
%
%      a(:,1) = lattice basis vector 1,
%      a(:,2) = lattice basis vector 2,
%      a(:,3) = lattice basis vector 3.
%
    
%% Outputs
%  vertex = 3XM array = M vertices of Wigner-Seitz cell
%  vertex(1:3,i) = Cartesian coordinates of vertex i
%
%  edge   = 2XN array = N edges of Wigner-Seitz cell
%  edge(1,n) = i <=> n-th edge starts at vertex(:,i)
%  edge(2,n) = j <=> n-th edge   ends at vertex(:,j)
%
%  face = NF*MV array where NF = number of faces 
%         of Wigner-Seitz cell and MV = maximum of the numbers of
%         vertices of each face and such that 
%
%         vertex(1:3,face(i,j)) = vertex i of face j
%
%  NB: if face(i,j) = NaN, this means that the number of
%      vertices of face j is smaller than i
%      (feature needed for plotting: see Matlab/Octave manual)
%  
%  vn = 5XNF array = lattice vectors normal to faces 
%       that are closest to the origin
%  vn(1:3,i) Cartesian coordinates of vector normal to face i
%  vn(4,i) = norm(v(1:3,i))
%  vn(5,i) = v(4,i)^2
%
 
%% Recognized options in varargin 
% (uppercases for readability are optional): 
    
% if varargin{k} = 'tol', 
%                  then varargin{k+1} = real scalar > 0 
%                  = tolerance to recognize zero.
%                  In practice, two floating-point numbers
%                  will be considered equal if the absolute value 
%                  of their difference is < tol   
    
% if varargin{k} = 'Origin', 
%                   then varargin{k+1} = [x,y,z]' 
%                   = coordinates of origin of Wigner-Seitz cell
    
% if varargin{k} = 'Space'
%                  then varargin{k+1} = 'Direct' or 'Recip'
%                  to work in direct (default) or reciprocal space.
%                  This option changes labels of axes and the color
%                  of the cell faces.
    
% if varargin{k} = 'FullPlot'
%                  then varargin{k+1} = Boolean 
%                  If true, the following plot options are set to
%                  true: 'PlotBasisVectors','PlotEdges','PlotFaces',
%                  'ShowNodesNormalToFaces'
    
% if varargin{k} = 'PlotBasisVectors'
%                  then varargin{k+1} = Boolean  
%                  If true, plot basis vectors of direct lattice  
    
% if varargin{k} = 'PlotEdges',
%                  then varargin{k+1} = Boolean controlling 
%                  If true, plot edges of Wigner-Seitz cell
%                  (polygons formed by edges are not filled
%                  with any color)
    
% if varargin{k} = 'PlotFaces',
%                  then varargin{k+1} = Boolean  
%                  If true, fill faces of Wigner-Seitz cell
%                  (polygons formed by edges are filled with
%                  colors on the basis of option 'ColorMap')
    
% if varargin{k} = 'Transparency', 
%                  then varargin{k+1} = real scalar > 0 
%                  = Transparency of faces 
%                    in range 0  (transparent) to 1 (opaque)

% if varargin{k} = 'PlotLatticeNodes', 
%                  then varargin{k+1} = Boolean  
%                  If true, plot of all lattice nodes
%                  used in the computation
    
% if varargin{k} = 'ShowNodesNormalToFaces', 
%                  then varargin{k+1} = Boolean 
%                  If true, show the lattice nodes normal to faces

% if varargin{k} = 'PlotCellVertices',
%                   then varargin{k+1} = Boolean  
%                   If true, plot Wigner-Seitz cell vertices 

% if varargin{k} = 'Verbose', 
%                  then varargin{k+1} = Boolean 
%                  If true, display some computation details  
    
% if varargin{k} = 'LatticeType', 
%                  then varargin{k+1 } = string
%                  =  Description of lattice type used for:
%                  [-] title of plot window;
%                  [-] switching fast procedure  
%                      for cubic lattices (3D only)
    
% if varargin{k} = 'ColorMap', 
%                  then varargin{k+1} = Character string 
%                  defining one of the pre-defined colormaps
%                  used to plot the faces of a cell
% List of MATLAB pre-defined colormaps:
% https://fr.mathworks.com/help/matlab/ref/colormap.html
% List of Octave pre-defined colormaps:
% https://octave.sourceforge.io/octave/function/colormap.html
    
% if varargin{k} = 'Colors', 
%                   then varargin{k+1 } = 7x3 matrix
%                   = Color order of successive objects
%                   (not including lattice nodes)
    
% if varargin{k} = 'ChangeFaceColor' (3D only), 
%                   then varargin{k+1} = Boolean
%                   If true, change color on each face 
    
% if varargin{k} = 'FaceColor' (3D only), 
%                   then varargin{k+1} = Integer in range 1 to 64
%                   = Fixed index in colormap for colouring faces
%                   This option overwrites the option 'ChangeFaceColor'

% if varargin{k} = 'ShowLegend', 
%                   then varargin{k+1} = Boolean
%                   If true, show legend   
    
% if varargin{k} = 'French', 
%                   then varargin{k+1} = Boolean
%                   If true, legends are in French 
    
%% INPUT PARSING
    
p=inputParser; p.CaseSensitive = false; % Initialize input parser

addRequired(p,'a', ...
             @(x) isreal(x) && all(size(x)==[3 3]));
addParameter(p,'tol',1e-12,...
             @(x) isreal(x) && isscalar(x) && x>0);
addParameter(p,'transparency',0.5,...
             @(x) isreal(x) && isscalar(x) && x>=0 && x<=1);
addParameter(p,'verbose',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'origin',[0 0 0]', ...
             @(x) isreal(x) && all(size(x)==[3 1]));
addParameter(p,'space','direct', ...
             @(x) ischar(x));
addParameter(p,'latticetype',' ', ...
             @(x) ischar(x));
addParameter(p,'fullplot',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plotfaces',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plotedges',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plotcellvertices',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plotbasisvectors',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'plotlatticenodes',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'shownodesnormaltofaces',false,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'showlegend',true,...
             @(x) islogical(x) && isscalar(x));
addParameter(p,'changefacecolor',false, ...
             @(x) islogical(x));
addParameter(p,'facecolor',53, ...
             @(x) isintegervalue(x) && x>0 && x<65);
addParameter(p,'colormap',[], ...
             @(x) ischar(x)); 
addParameter(p,'colors',get(gca,'colororder'),...
             @(x) isreal(x) && all(size(x)==[7 3]));
addParameter(p,'french',false,...
             @(x) islogical(x) && isscalar(x));

parse(p,a,varargin{:}); % Read the argument list knowing the rules
a=p.Results.a;
origin=p.Results.origin;
space=p.Results.space;
verbose=p.Results.verbose;
tol=p.Results.tol;
ty=p.Results.transparency;
lattice_type=p.Results.latticetype;
full_plot=p.Results.fullplot;
plot_faces=p.Results.plotfaces;
plot_edges=p.Results.plotedges;
plot_cell_vertices=p.Results.plotcellvertices;
plot_basis_vectors=p.Results.plotbasisvectors;
plot_lattice_nodes=p.Results.plotlatticenodes;
show_nodes_normal_to_faces=p.Results.shownodesnormaltofaces;
show_legend=p.Results.showlegend;
changefacecolor=p.Results.changefacecolor;
facecolor=p.Results.facecolor;
color_map=p.Results.colormap; 
color_scheme=p.Results.colors;
French=p.Results.french;

if full_plot
    plot_faces=true; plot_edges=true;
    plot_basis_vectors=true;
    show_nodes_normal_to_faces=true;
end

%% LABELS OF WINDOW, AXES AND OF BASIS VECTORS

switch lower(space)
  case {'direct'}
    Title='Direct space: Wigner-Seitz cell';
    coord{1}='\itx\rm_1';coord{2}='\itx\rm_2';coord{3}='\itx\rm_3';
    bvec{1} ='\bfa\rm_1'; bvec{2}='\bfa\rm_2'; bvec{3}='\bfa\rm_3';
    bvec{4}='\bfa\rm\it_j';
    if isempty(color_map)
        color_map='viridis'; facecolor=32;
    end
  case {'recip','reciprocal'}
    Title='Reciprocal space: Brillouin Zone (BZ)';
    coord{1}='\itq\rm_1';coord{2}='\itq\rm_2';coord{3}='\itq\rm_3';
    bvec{1} ='\bfg\rm_1'; bvec{2}='\bfg\rm_2'; bvec{3}='\bfg\rm_3';
    bvec{4}='\bfg\rm\it_j';
    if isempty(color_map)
        color_map='pink'; facecolor=53;
    end
  otherwise
    error('function WignerSeitz: unrecognized parameter space=%s',space)
end

%% DETERMINE DIMENSION

for i=1:3
    norma(i)=norm(a(:,i));
end
dimension=sum(norma>tol);

if verbose
    fprintf('\n*************************************************\n')
    fprintf(' WIGNER-SEITZ CELL OF ORDER %d IN %dD %s SPACE\n',...
            1,dimension,upper(space))
    fprintf('*************************************************\n')
end

%% 1-DIMENSION

if (dimension==1) 
    
    [R,neighbour]=lattice(a,2);
    
    k=0; 
    for j=1:3
        if (norma(j) > tol)     % 1-DIMENSION ALONG a(:,j)
            i=j; 
        else
            k=k+1; vn(1:3,k)=a(1:3,j);   % Vectors normal to edge
            vn(5,k)=v(1,k)^2+v(2,k)^2+v(3,k)^2; v(4,1)=sqrt(v(5,k));
        end
    end 
    vertex(:,1)= -a(:,i)/2;
    vertex(:,2)=  a(:,i)/2;
    edge(1,1)=1; edge(2,1)=2; 
    fprintf('Function WignerSeitz in %s space: 2 vertices, 1 edge.\n',...
            space);
    face=[1 2 1];  % Face reduced to an edge    
end

%% 2-DIMENSION

if (dimension==2) 
    if (norm(a(:,3)) < tol )     % 2-DIMENSION IN PLANE x_1-x_2
        ii=1; jj=2; iz=3;
    end
    if (norm(a(:,1)) < tol )     % 2-DIMENSION IN PLANE x_2-x_3
        ii=2; jj=3; iz=1;
    end
    if (norm(a(:,2)) < tol )     % 2-DIMENSION IN PLANE x_1-x_3
        ii=1; jj=3; iz=2;
    end
    
    %% VECTORS NORMAL TO FACE
    
    TU=eye(3,3);
    vn(1:3,1)=TU(1:3,iz); vn(1:3,2)=-TU(1:3,iz); 
    for k=1:2
        vn(5,k)=vn(1,k)^2+vn(2,k)^2+vn(3,k)^2; vnn(4,1)=sqrt(vn(5,k));
    end
    
    %% WORK SET OF 2D LATTICE VECTORS
    
    [T,R,~]=lattice(a,1,1); 
    
    %% Projection to plane x_i-x_j
    
    for n=1:columns(T)
        TT(1,n)=T(ii,n); 
        TT(2,n)=T(jj,n);
    end
    
    for n=1:columns(R)
        RR(1,n)=R(ii,n); 
        RR(2,n)=R(jj,n);
    end
    
    %% Optimization of the work set RR
    %  
    %  If a vector TT(:,m) is not colinear 
    %  with any of the vectors of the current work set RR, 
    %  this means that the direction TT(:,m) is unexplored by the
    %  work set => update work set RR by including TT(:,m).
    %
    %  If a vector TT(:,m) is colinear with some vector RR(:,n)
    %  of the current work set RR and pointing in the same
    %  direction as RR(:,n), there is no need :
    %  [-] to update the work set RR;
    %  [-] to consider the TT(:,j) that are beyond (as viewed from
    %      the origin) the line normal to RR(:,n) and containing RR(:,n).
    
    nodes=columns(RR); newmax=columns(RR); exclus=[];
    for m=nodes+1:columns(TT)    % Search beyond first neighbours
        if (~ismember(m,exclus)) % Check if excluded by a previous step
            colinear=false;
            for n=1:newmax
                sp=TT(1,m)*RR(1,n)+TT(2,m)*RR(2,n);
                cosinus=sp/( sqrt(TT(1,m)^2+TT(2,m)^2) * sqrt(RR(1,n)^2+RR(2,n)^2) );
                if (abs(cosinus-1)<tol)
                    colinear=true; % TT(:,m) and RR(:n) are colinear and
                    break;         % pointing in the same direction
                end
            end
 
            if (colinear)
                for j=m+1:columns(TT)
                    % prp = position of TT(:,j) relatively to plane normal
                    % to RR(:,n) and including RR(:,n)
                    prp=R(1,n)*T(1,j)+R(2,n)*T(2,j)-RR(1,n)^2-RR(2,n)^2;
                    if (prp>0) % TT(:,j) beyond plane (as viewed from the origin)
                        exclus = [exclus j]; % Excluding TT(;,j)
                    end
                end
            else
                newmax=columns(RR)+1; % If TT(:,m) not in the same direction
                RR(:,newmax)=TT(:,m); % as any vector of the current work
                                      % set, include TT(:,m) in work set
            end
        end
    end
   
    fprintf(['Function WignerSeitz using %d nodes after searching ' ...
             'colinear vectors from the origin.\n'],columns(RR))
        
    %% VERTICES OF THE 2D WIGNER-SEITZ CELL
    
    c=RR/2; % The lines of the Wigner-Seitz unit cell pass through
                 % a subset of the points c=zone*RR/2
    
    for i=1:columns(RR) 
        d(i) = RR(1,i)*c(1,i) + RR(2,i)*c(2,i);
    end
    
    % The equation of a line orthogonal to RR(1:2,i) and containing
    % c(1:2,i) = RR(1:2,i)/2 reads:
    %
    % RR(1,i)*x(1) + RR(2,i)*x(2) = d(i) 
    % where
    % d(i) = RR(1,i)*c(1,i) + RR(2,i)*c(2,i)
    
    % The intersection between 2 lines i,j 
    % = a unit cell vertex = solution of the 2X2 linear system of equations:
    %
    % RR(1,i)*x(1) + RR(2,i)*x(2) = d(i) 
    % RR(1,j)*x(1) + RR(2,j)*x(2) = d(j) 
     
    
    %% Search possible vertices
    
    n=0;
    for i=2:columns(RR) % Consider all possible intersections of two lines
        for j=i+1:columns(RR)  
            b = [d(i), d(j)]';
            S(1,1:2) = RR(1:2,i);
            S(2,1:2) = RR(1:2,j);
            if( abs(det(S)) > tol ) 
                x = S\b; % Solve 2X2 linear system of equations
                n=n+1; u(1:2,n)=x(1:2);
            end
        end
    end
    
    if (verbose) 
        fprintf(['Function WignerSeitz: %d intersections as possible ' ...
                 'vertices.\n'],columns(u));
    end
    
    %% Search intersections closest to the origin 
    %% than to any other lattice vector

    k=0; 
    for i=1:n
        dorigin = u(1,i)^2 + u(2,i)^2 ; flag=true;
        for j=2:columns(RR)
            dr = (u(1,i)-RR(1,j))^2 + (u(2,i)-RR(2,j))^2;
            if( dr < dorigin-tol)
                flag=false; % Discard intersection closer to any lattice
                            % vector than to the origin
            end
        end
        if flag 
            k=k+1; v(1:2,k)=u(1:2,i);  
        end
    end
    
    if (verbose)  
        fprintf(['Function WignerSeitz: %d intersections closest ' ...
                 'to origin than to any other lattice vector.\n'],columns(v));
    end
    
    %% Discard eventual vertex duplicates
    
    k=0;
    for i=1:columns(v)
        check=true;
        for j=i+1:columns(v)
            if ( abs(v(1,i)-v(1,j)) < tol && abs(v(2,i)-v(2,j)) < tol)
                check=false; % Discard duplicate
            end
        end
        if (check) 
            k=k+1; vertex(1:2,k)=v(1:2,i); vertex(3,k)=0;
        end
    end
    
    if verbose 
        fprintf(['Function WignerSeitz: %d vertices after discarding ' ...
                 'duplicates.\n'], columns(vertex));
    end
    
    %% BACK TO 3-DIMENSION VECTORS
      
    buf=zeros(3,columns(vertex));
    for n=1:columns(vertex)
        buf(ii,n)=vertex(1,n); buf(jj,n)=vertex(2,n);
    end
    vertex = buf;
    
    for n=1:columns(RR)
        R(ii,n)=RR(1,n); R(jj,n)=RR(2,n); R(iz,n)=0;
    end
   
    %% EDGES OF 2D WIGNER-SEITZ CELL
     
    ptxy=vertex; ptxy(iz,:)=[];
    ppt=Jarvis(ptxy); 
    % Function Jarvis returns the convex hull of
    % a set of points in a plane, thereby providing the correct 
    % rotation order of the vertices defining a face 
    % that will be suitable for plotting by the function
    % "patch" (see Matlab/Octave documentation).
    
    for n=1:columns(vertex)
        vertex(ii,n)=ppt(1,n); vertex(jj,n)=ppt(2,n);
    end
    
    for n=1:columns(vertex) % Number of edges of a polygon
                            % = number of its vertices
        edge(1,n)=n; edge(2,n)=n+1;
        if (n==columns(vertex))
            edge(2,n)=1;
        end
        
    end
    
    %% Vectors normal to edges
    
    for n=1:columns(vertex)-1
        vne(1:3,n)=vertex(1:3,n)+vertex(1:3,n+1);
    end
    vne(1:3,columns(vertex))=vertex(1:3,columns(vertex))+vertex(1:3,1);
    
    %% SINGLE FACE OF 2D WIGNER-SEITZ CELL
    
    face=[ 1:1:columns(vertex) ];
    
    fprintf('Function WignerSeitz in 2-D %s space: %d vertices, %d edges.\n',...
            space,columns(vertex), columns(vertex));

end % End of 2-dimension
    
%% 3-DIMENSION

if (dimension==3)
    
    %% WORK SET OF 3D LATTICE VECTORS
    
    switch lattice_type
      case {'cP','cI','cF'}
        [~,R,~]=lattice(a,1,2); % Fast lane for high symmetry lattices
        T=R;
      otherwise
        [T,R,~]=lattice(a,3,2); % R(:,n) is the starting work set that
                                % will expanded if needed using 
                                % the "provision" of vectors in array T
    end
    
    %% Optimization of the work set R
    %  
    %  If a vector T(:,m) is not colinear 
    %  with any of the vectors of the current work set R, 
    %  this means that the direction T(:,m) is unexplored by the
    %  work set => update work set R by including T(:,m).
    %
    %  If a vector T(:,m) is colinear with some vector R(:,n)
    %  of the current work set R and pointing in the same
    %  direction as R(:,n), there is no need :
    %  [-] to update the work set R;
    %  [-] to consider the T(:,j) that are beyond (as viewed from
    %      the origin) the plane normal to R(:,n) and containing R(:,n).
    
    nodes=columns(R); newmax=columns(R); exclus=[];
    for m=nodes+1:columns(T)     % Search beyond first set of nodes
        if (~ismember(m,exclus)) % Check if excluded by a previous step
            colinear=false;
            for n=1:newmax
                sp=T(1,m)*R(1,n)+T(2,m)*R(2,n)+T(3,m)*R(3,n);
                cosinus=sp/(T(4,m)*R(4,n));
                if (abs(cosinus-1)<tol)
                    colinear=true; % T(:,m) and R(:n) are colinear and
                    break;         % pointing in the same direction
                end
            end
 
            if (colinear)
                for j=m+1:columns(T)
                    % prp = position of T(:,j) relatively to plane normal
                    % to R(:,n) and including R(:,n)
                    prp=R(1,n)*T(1,j)+R(2,n)*T(2,j)+R(3,n)*T(3,j)-R(5,n);
                    if (prp>0) % T(:,j) beyond plane (as viewed from the origin)
                        exclus = [exclus j]; % Excluding T(;,j)
                    end
                end
            else
                newmax=columns(R)+1;% If T(:,m) not in the same direction
                R(:,newmax)=T(:,m); % as any vector of the current work
                                    % set, include T(:,m) in work set
            end
        end
    end
    if(verbose)
        fprintf(['Function WignerSeitz using %d nodes after searching ' ...
                 'colinear vectors from the origin.\n'],columns(R))
    end
    
    %% VERTICES OF 3D WIGNER-SEITZ CELL
    
    c=R/2; % The planes of the Wigner-Seitz unit cell pass through
           % a subset of the points c=R/2
    
    for i=1:columns(R) 
        d(i) = R(1,i)*c(1,i) + R(2,i)*c(2,i)+ R(3,i)*c(3,i);
    end
    
    % The equation of a plane orthogonal to R(1:3,i) and containing
    % c(1:3,i) = R(1:3,i)/2 reads:
    %
    % R(1,i)*x(1) + R(2,i)*x(2) + R(3,i)*x(3) = d(i) 
    % where
    % d(i) = R(1,i)*c(1,i) + R(2,i)*c(2,i)+ R(3,i)*c(3,i)
    
    % The intersection between 3 planes i,j,k 
    % = a unit cell vertex = solution of the 3X3 linear system of equations:
    %
    % R(1,i)*x(1) + R(2,i)*x(2) + R(3,i)*x(3) = d(i) 
    % R(1,j)*x(1) + R(2,j)*x(2) + R(3,j)*x(3) = d(j) 
    % R(1,k)*x(1) + R(2,k)*x(2) + R(3,k)*x(3) = d(k) 
    
    %% Search possible vertices
    
    nodes=columns(R); nmax=(nodes-1)*(nodes-2)*(nodes-3);
    u=zeros(3,nmax);    
    
    n=0; 
    if(verbose)
        fprintf('Function WignerSeitz searching intersections ...');
    end
    for i=2:nodes
        for j=i+1:nodes 
            for k=j+1:nodes
                b = [d(i), d(j), d(k)]';
                S(1,1:3) = R(1:3,i);
                S(2,1:3) = R(1:3,j);
                S(3,1:3) = R(1:3,k);
                if( abs(det(S)) > tol ) 
                    x = S\b; % Solve 3X3 linear system of equations
                    n=n+1; u(1:3,n)=x(1:3);
                end
            end
        end
    end
    u(:,n+1:nmax)=[];
    if (verbose) 
        fprintf([' %d intersections as possible ' ...
                 'vertices.\n'],n);
    end
    
    %% Search intersections closest to the origin 
    %% than to any other lattice vector
    
    k=0; v=zeros(3,columns(u));
    if(verbose)
        fprintf(['Function WignerSeitz searching intersections closest ' ...
                 'to origin...\n']);  
    end
    
    for i=1:n
        flag=true; dorigin = u(1,i)^2 + u(2,i)^2 + u(3,i)^2;
        for j=2:columns(R)
            dr = (u(1,i)-R(1,j))^2 + (u(2,i)-R(2,j))^2 + (u(3,i)-R(3,j))^2;
            if (dr < dorigin-tol)
                flag=false; % Discard intersection closer to
                break;      % any lattice vector than to the origin
            end
        end
        if (flag)
            k=k+1; v(1:3,k)=u(1:3,i);  
        end
    end
    v(:,k+1:n)=[];
    if verbose  
        fprintf(['Function WignerSeitz: %d intersections closest ' ...
                 'to origin than to any other lattice vector.\n'],columns(v));
    end
    
    %% Discard eventual vertex duplicates
    
    k=0;
    for i=1:columns(v)
        check=true;
        for j=i+1:columns(v)
            if ( abs(v(1,i)-v(1,j)) < tol && abs(v(2,i)-v(2,j)) < tol && ...
                 abs(v(3,i)-v(3,j)) < tol)
                check=false; % Discard duplicate
            end
        end
        if (check) 
            k=k+1; vertex(1:3,k)=v(1:3,i); 
        end
    end
    if (verbose) 
        fprintf(['Function WignerSeitz: %d vertices after discarding ' ...
                 'duplicates.\n'], columns(vertex));
    end
    
    %% EDGES OF WIGNER-SEITZ CELL
    
    % The edges of the unit cell are found by considering all
    % pairs of vertices. If a pair of vertices share two planes of the
    % Wigner-Seitz cell boundaries, there is an edge between these two vertices.
    
    n=0;
    for i=1:columns(vertex)
        for j=i+1:columns(vertex)
            for k=2:columns(R) 
                % w(1) = 0 if vertex i in plane k
                w(1) = R(1,k)*vertex(1,i)+R(2,k)*vertex(2,i)+...
                       R(3,k)*vertex(3,i)-d(k);
                % w(2) = 0 if vertex j in plane k
                w(2) = R(1,k)*vertex(1,j)+R(2,k)*vertex(2,j)+...
                       R(3,k)*vertex(3,j)-d(k); 
                if (abs(w(1))<tol && abs(w(2))<tol) 
                    for l=k+1:columns(R) 
                        % w(3) = 0 if vertex i in plane l
                        w(3) = R(1,l)*vertex(1,i)+R(2,l)*vertex(2,i)+...
                               R(3,l)*vertex(3,i)-d(l); 
                        % w(4) = 0 if vertex j in plane l
                        w(4) = R(1,l)*vertex(1,j)+R(2,l)*vertex(2,j)+...
                               R(3,l)*vertex(3,j)-d(l);
                        if(abs(w(3))<tol && abs(w(4))<tol)% Edge connecting
                            n=n+1;edg(1,n)=i;edg(2,n)=j;  % vertices i & j     
                        end
                    end
                end
            end
        end   
    end
    
    %% Discard eventual edge duplicates
    
    k=0;
    for i=1:columns(edg)
        check=true;
        for j=i+1:columns(edg)
            if ( abs(edg(1,i)-edg(1,j)) < tol && abs(edg(2,i)-edg(2,j)) < tol)
                check=false; % Discard duplicate
            end
        end
        if check 
            k=k+1; edge(1:2,k)=edg(1:2,i); 
        end
    end
    
    if verbose
        fprintf('Function WignerSeitz: %d edges have been found.\n',...       
                 columns(edge))
    end
        
    %% FACES OF WIGNER-SEITZ CELL
    
    % The number of faces is known by Euler's characteristic 
    % formula for convex polyhedrons
    
    Euler=2+columns(edge)-columns(vertex);
    
    %% The faces of the unit cell are found
    %% if a plane contains more than two vertices.

    nmax=0; kf=0; face=zeros(columns(R),columns(vertex));
    for k=2:columns(R) % Avoid starting with R=[0,0,0]
        for i=1:columns(vertex)
            n=0; buf=[];
            for j=1:columns(vertex)
                % w(1) = 0 if vertex i in plane k
                w(1) = R(1,k)*vertex(1,i)+R(2,k)*vertex(2,i)+...
                       R(3,k)*vertex(3,i)-d(k);
                % w(2) = 0 if vertex j in plane k
                w(2) = R(1,k)*vertex(1,j)+R(2,k)*vertex(2,j)+...
                       R(3,k)*vertex(3,j)-d(k);            
                if ( abs(w(1))<tol && abs(w(2))<tol )
                    n=n+1; buf(n)=j;  
                end
            end
            if(n>2) % More than 2 vertices in plane k 
                    % Add to face list if not already found
                if (kf==0 || (kf>0 && ~ismember(k,vnind))) 
                    kf=kf+1;
                    face(kf,1:n)=buf(1:n);
                    vnind(kf)=k; % Store index of vector normal to face
                    if (n>nmax)
                        nmax=n;
                    end
                end
            end
        end
    end
      
    face(:,nmax+1:end)=[]; face(kf+1:end,:)=[]; % Minimize array shape 
    
    if(rows(face) ~= Euler)
        fprintf(['\nWARNING: inconsistency in function WignerSeitz: number of ' ...
                 'faces = %d not equal to Euler characteristics = %d\n'],...
                rows(face),Euler)
    else
        fprintf(['Function WignerSeitz in 3D %s space: %d vertices, ' ...
                 '%d edges, %d faces.\n'],...
                 space,columns(vertex),columns(edge),Euler);
    end

    for i=1:rows(face)
        kept(i)=columns(face);
        for j=columns(face):-1:1
            if(face(i,j)==0)
                face(i,j) = NaN; kept(i)=j-1; 
                % See Matlab/Octave manual: NaN allows to plot faces with
                % different numbers of vertices
            end
        end
    end    
 
    %% ORDERING VERTICES OF EACH FACE AS A CONVEX POLYGON FOR PLOTTING
    
    for i=1:rows(face) % Number of faces = number of vectors normal
                       % to faces
        vn(1:5,i)=R(1:5,vnind(i));    % Lattice vector normal to face
        theta= acos(vn(3,i)/vn(4,i)); % In spherical coordinates, angles
        phi  = atan2(vn(2,i),vn(1,i));% of vector normal to face i
                                        
        ROTZ=RotationAxe(3,-phi);   % Rotation around axe 3
        ROTY=RotationAxe(2,-theta); % Rotation around axe 2
        
        ROTATION=ROTY*ROTZ; % Rotation bringing vn(:,i) parallel to axe 3
       
        pt=[]; ptxy=[]; ppt=[];
        for j=1:kept(i)     % Rotation  applied to all vertices
            pt(1:3,j)=ROTATION*[vertex(1:3,face(i,j))];
        end
        % The resulting pt(:,j) have all the same third component
        % (the rotation moved them in a plane parallel to the plane
        % of Cartesian basis vectors 1 and 2).
        % This third component must be stored to allow a correct 
        % inverse rotation later. Therefore pt is duplicated 
        % as ptxy of which we keep only the components 
        % in the plane of axes 1 and 2 as arguments of the 
        % Jarvis' march algorithm.
        ptxy=pt; ptxy(3,:)=[];
       
        ppt=Jarvis(ptxy); 
        % Function Jarvis returns the convex hull of
        % a set of points in a plane, thereby providing the correct 
        % rotation order of the vertices defining a face 
        % that will be suitable for plotting by the function
        % "patch" (see Matlab/Octave documentation).
        
        for j=1:kept(i)
            pt(1,j)=ppt(1,j); pt(2,j)=ppt(2,j); 
            % pt is now in the correct order to define a convex
            % polygon for the current face. The third component
            % was not changed because it is identical for
            % all points after the rotation of axes.
        end
        
        ROTZ=RotationAxe(3,phi);   % Inverse rotation around axe 3
        ROTY=RotationAxe(2,theta); % Inverse rotation around axe 2
        
        ROTATION=ROTZ*ROTY; % Inverse rotation bringing axe 3 back in the
                            % direction of vn(:,i) 
        
        for j=1:kept(i) % Loop over each vertex of face i
            pto=ROTATION*[pt(1:3,j)]; % Re-ordered vertex after inverse rotation
            for k=1:columns(vertex)   % Search match of pto in array vertex(:,k)
                vk=[vertex(1:3,k)];   % Buffer
                if (norm(pto-vk)< tol)% Matching...
                    face(i,j)=k;      % ... modifies vertex j of face i
                end
            end
        end
    end  
    if verbose
        fprintf(['Function WignerSeitz ordered vertices for plotting ' ...
                 'faces by function patch.\n']) 
        fprintf('\n***** End of function WignerSeitz ***************\n')
    end
end % End of 3-dimension
 

%% ANYTHING TO PLOT ?

if ~any([plot_faces ...
         plot_edges ...
         plot_cell_vertices ...
         plot_basis_vectors ...
         plot_lattice_nodes ...
         show_nodes_normal_to_faces])
    return;
end

%% OPTIONAL PLOTS
    
plotopt=plot_preamble; 
k=plotopt.k; Legend=plotopt.Legend;
color_scheme=plotopt.color_scheme;
lw=plotopt.lw;   elw=plotopt.elw;
mhs=plotopt.mhs; mks=plotopt.mks;

if (plot_basis_vectors) 
    for i=1:3
        k=k+1; Legend{k}=bvec{i};
        quiver3(origin(1),origin(2),origin(3),a(1,i),a(2,i),a(3,i),...
                'color',color_scheme(i,1:3),...
                'linewidth',lw,...
                'maxheadsize',mhs,...
                'DisplayName',Legend{k}); hold on;
    end
end

if (plot_lattice_nodes) 
    k=k+1; 
    if French
        Legend{k}='Noeuds du réseau';
    else
        Legend{k}='Lattice nodes';
    end
    plot3(R(1,:),R(2,:),R(3,:),'o',...
          'MarkerSize',mks,...
          'Color',color_scheme(3,:),...
          'DisplayName',Legend{k}); hold on;
end


if (show_nodes_normal_to_faces)
    if (dimension == 2)
        k=k+1; 
        if French
            Legend{k}='Noeuds orthogonaux aux arêtes';
        else 
            Legend{k}='Lattice nodes normal to edges';
        end
        plot3(vne(1,:),vne(2,:),vne(3,:),'o',...
              'MarkerSize',mks,...
              'MarkerFaceColor',color_scheme(2,:),...
              'Color',color_scheme(3,:),...
              'DisplayName',Legend{k}); hold on;
    end
    if (dimension == 3)
        k=k+1; 
        if French
            Legend{k}='Noeuds orthogonaux aux faces'; 
        else
            Legend{k}='Lattice nodes normal to faces';
        end
        plot3(vn(1,:),vn(2,:),vn(3,:),'o',...
              'MarkerSize',mks,...
              'Color',color_scheme(3,:),...
              'MarkerFaceColor',color_scheme(2,:),...
              'DisplayName',Legend{k}); hold on;
    end
end

if (plot_cell_vertices) 
    k=k+1; 
    if French
        Legend{k}='Sommets de la maille de Wigner-Seitz';
    else
        Legend{k}='Vertices of Wigner-Seitz cell';
    end
    plot3(vertex(1,:),vertex(2,:),vertex(3,:),'o',...
          'color',color_scheme(2,:),...
          'DisplayName',Legend{k}); hold on;
end

if (plot_edges)  
    if (~plot_faces)
        k=k+1;
        if French
            Legend{k}='Arêtes de la maille de Wigner-Seitz'; 
        else
            Legend{k}='Edges of Wigner-Seitz cell'; 
        end
        pts = [vertex(:,edge(1,1))'; vertex(:,edge(2,1))'];
        plot3(pts(:,1), pts(:,2), pts(:,3),...
              '-','color',color_scheme(7,:),...
              'DisplayName',Legend{k}); hold on; 
        for i=2:columns(edge)
            pts = [vertex(:,edge(1,i))'; vertex(:,edge(2,i))'];
            plot3(pts(:,1), pts(:,2), pts(:,3),...
                  '-','color',color_scheme(7,:),...
                  'HandleVisibility','off'); hold on; 
        end
    end
end

if (plot_faces) 
    k=k+1; 
    switch lower(space)
      case {'direct'}
        if French
            Legend{k}='Maille de Wigner-Seitz'; 
        else
            Legend{k}='Wigner-Seitz cell'; 
        end
      case {'recip','reciprocal'}
        if French
            Legend{k}='Zone de Brillouin'; 
        else
            Legend{k}='Brillouin zone'; 
        end
    end
    
    rgb=colormap(color_map);
    
    % Set colormap size to 64
    % (Octave's default, whereas Matlab's default = 256)
    P = size(rgb,1); N = min(64,P);
    cmap= interp1(1:P, rgb, linspace(1,P,N), 'linear');
    
    switch dimension % Color and transparency of faces
      case 1          
        color_face(1,1:3)=color_scheme(4,:);      
      case 2
        color_face(1,1:3)=color_scheme(4,:);
      case 3
        % Try to avoid extreme colors when changing face.
        % Recommended colormaps:
        % viridis,gray,bone,pink,autumn,ocean,winter
        if changefacecolor
            for i=1:rows(face)
                color_face(i,1:3)=cmap(53-i*4,1:3);
            end
        else % Fixed color
            for i=1:rows(face)
                color_face(i,1:3)=cmap(facecolor,1:3);
            end
        end
    end
    
    patch('Vertices', vertex', 'Faces', face, ...
          'LineStyle','-',...
          'LineWidth', elw,...
          'EdgeColor',color_scheme(7,:),... % Options: 'none', 'flat', 'interp', [R G B]
          'EdgeAlpha',1,...     % Edge transparency in range 0 (transparent) to 1 (opaque)
          'FaceColor','flat',...% Options: 'none', 'flat', 'interp', [R G B]
          'FaceAlpha',ty,...    % Face transparency in range 0  (transparent) to 1 (opaque)
          'FaceVertexCData',color_face,...
          'DisplayName',Legend{k}); hold on;
    
end

%% AXES, LABELS & LEGEND

axis equal; grid on; legend(Legend,'Location','EastOutside');
xlabel(coord{1}); ylabel(coord{2}); zlabel(coord{3});     

if(~show_legend)
    legend(gca,'off');
end

%% POINT OF VIEW 

switch dimension
  case 1
    pv=[0,0,1];  
    if (norm(a(:,3)) > tol)
      pv=[0,1,0];
    end
    
  case 2
    pv=[0,0,0];  
    for k=1:3
        if (norm(a(:,k)) < tol)
            pv(k)=-((-1)^k);
        end
    end
    
  case 3
    pv=[1,0.5,0.3];
 
end
view(pv);

end % End of function WignerSeitz
    