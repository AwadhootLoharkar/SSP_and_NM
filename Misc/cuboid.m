function cuboid(a,varargin)

%%  3-DIMENSIONAL PLOT OF THE CUBOID HEXAHEDRON
%%  SPANNED BY 3 BASIS VECTORS

%%  Input [length]:

%  a = 3X3 array of basis vectors in column order
%      i.e: a(1:3,i) = basis vector i 
%      in any unit of length.
    
%% Output
    
%  3-dimensional plot of cuboid hexahedron spanned by basis vectors 
    
%% Recognized options in varargin 
% (uppercases for readability are optional): 

% if varargin{k} = 'ColorMap', 
%                   then varargin{k+1} = Character string 
%                   defining one of the pre-defined colormaps
    
% if varargin{k} = 'ChangeFaceColor'  
%                   then varargin{k+1} = Boolean
%                   If true, change color on each face 
    
% if varargin{k} = 'FaceColor' (3D only), 
%                   then varargin{k+1} = Integer in range 1 to 64
%                   = Fixed index in colormap for colouring faces
%                   This option overwrites the option 'ChangeFaceColor'

% if varargin{k} = 'Transparency', 
%                  then varargin{k+1} = real scalar > 0 
%                  = Transparency of faces 
%                    in range 0  (transparent) to 1 (opaque)
    
% if varargin{k} = 'DisplayName', 
%                   then varargin{k+1} = Character string 
%                   that will appear in legend box
    
% if varargin{k} = 'Origin', 
%                   then varargin{k+1} = [x,y,z]' 
%                   = coordinates of origin of the basis vectors
    
% if varargin{k} = 'ShowLegend', 
%                   then varargin{k+1} = Boolean 
%                   If true, show legend  

% List of MATLAB pre-defined colormaps:
% https://fr.mathworks.com/help/matlab/ref/colormap.html
    
% List of Octave pre-defined colormaps:
% https://octave.sourceforge.io/octave/function/colormap.html

%% INPUT PARSING RULES 
p=inputParser; p.CaseSensitive = false;
addRequired(p,'a', @(x) isreal(x) && all(size(x)==[3 3])) ;  
addParameter(p,'origin',[0 0 0]',...
             @(x) isreal(x) && all(size(x)==[3 1]));
addParameter(p,'transparency',0.2,...
             @(x) isreal(x) && isscalar(x) && x>=0 && x<=1);
addParameter(p,'colormap','viridis', ...
             @(x) ischar(x)); 
addParameter(p,'changefacecolor',false, ...
             @(x) islogical(x));
addParameter(p,'facecolor',42, ...
             @(x) isintegervalue(x) && x>0 && x<65);
addParameter(p,'showlegend',true,...
             @(x) islogical(x) && isscalar(x)); 
addParameter(p,'displayname','Cuboid spanned by basis vectors',...
             @(x) ischar(x)); 

%% PARSE INPUT
parse(p,a,varargin{:}); % Read the argument list knowing the rules
a=p.Results.a;           
origin=p.Results.origin;
ty=p.Results.transparency;
color_map=p.Results.colormap;
changefacecolor=p.Results.changefacecolor;
facecolor=p.Results.facecolor;
show_legend=p.Results.showlegend;
display_name=p.Results.displayname;

%% CORE JOB

vertex(1:3,1) = origin(1:3);
vertex(1:3,2) = origin(1:3) + a(1:3,1);
vertex(1:3,3) = origin(1:3) + a(1:3,1) + a(1:3,2);
vertex(1:3,4) = origin(1:3) + a(1:3,2);
vertex(1:3,5) = origin(1:3) + a(1:3,3);
vertex(1:3,6) = origin(1:3) + a(1:3,1) + a(1:3,3);
vertex(1:3,7) = origin(1:3) + a(1:3,1) + a(1:3,2) + a(1:3,3);
vertex(1:3,8) = origin(1:3) + a(1:3,2) + a(1:3,3);

face = [ 1, 2, 3, 4;
         5, 6, 7, 8;
         1, 2, 6, 5;
         4, 3, 7, 8;
         2, 3, 7, 6;
         4, 1, 5, 8 ];

%% PLOT

plotopt=plot_preamble; k=plotopt.k; Legend=plotopt.Legend;

k=k+1; Legend{k}=display_name;

rgb=colormap(color_map);

% Set colormap size to 64
% (Octave'default, whereas Matlab's default = 256)
P = size(rgb,1); N = min(64,P);
cmap= interp1(1:P, rgb, linspace(1,P,N), 'linear');

% Try to avoid extreme colors when changing face.
% Recommended colormaps:
% viridis,gray,bone,pink,autumn,ocean,winter
if changefacecolor
    for i=1:6
        color_face(i,1:3)=cmap(53-i*4,1:3);
    end
else
    for i=1:6
        color_face(i,1:3)=cmap(facecolor,1:3);
    end
end

patch ('Vertices', vertex', 'Faces', face, ...
       'LineStyle','-',...
       'LineWidth', plotopt.elw,...
       'EdgeColor',[0 0 0],...% Options: 'none', 'flat', 'interp', [R G B]
       'EdgeAlpha',1,...      % Edge transparency in range 0 (transparent) to 1 (opaque)
       'FaceColor', 'flat',...% Options: 'none', 'flat', 'interp', [R G B]
       'FaceAlpha',ty,...    % Face transparency in range 0 (transparent) to 1 (opaque)
       'FaceVertexCData', color_face,...
       'DisplayName',Legend{k}); hold on;

% Axes, Labels & Legend
    
axis equal; grid on; legend(Legend,'Location','EastOutside');

if(~show_legend)
    legend(gca,'off');
end

end % End of function cuboid
    