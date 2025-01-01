function [R]=digitizer(inputfile,xcalib,ycalib,varargin)

%% digitizer.m : DIGITIZE A CURVE IN AN IMAGE FILE
%
%  External functions not built in Matlab/Octave: 
%  changefiletype, printtable
%
%% The three first points MUST be calibration points
%  of the graph, namely:
%
%  point 1 = origin of the graph to be digitized
%  point 2 = the abscissa tick mark featuring xcalib
%            as abscissa in graph unit
%  point 3 = the ordinate tick mark featuring ycalib
%            as ordinate in graph unit
%
%% After calibration:
%
%  points >3 Data points that may be entered randomly.
%            They will be sorted in order of increasing abscissae. 

%  NB: On the image, all points will be marked by a green spot.
%      The spots marking calibration points will be inside blue
%      squares.
%      The spots marking data points will be inside red squares.

%% ARGUMENTS
%
%  inputfile = Image file name structured as fname.ftype
%              where ftype is an image format 
%              that can be handled by imread, 
%              e.g.: gif, jpeg, png, pbm, etc... BUT NOT pdf!
%
%  xcalib = abscissa of expected point 2 in graph unit
%  ycalib = ordinate of expected point 3 in graph unit

%% MOUSE INPUT
%
%  Button 1: Add point; Button 2: Remove last point;
%  Button 3: Close image = end of data acquisition

%% OUTPUT TO FILE
%
%  R = 2XN array of data points in graph units (stored column-wise)
%      (excluding calibration points) and such that (with n=1:N):
%      R(1,n) = abscissa of point n  
%      R(2,n) = ordinate of point n  

%% OUTPUT FILENAME
%  = <inputfile> whose extension is changed to .dat (default)
%  = <outputfile> given as argument of option 'OutputFile'
%    (see below: Recognized options in varargin)

%% Recognized options in varargin 
% (uppercases for readability are optional): 
%
% if varargin{k} = 'OutputFile', 
%                  then varargin{k+1} = Character string 
%                  = name of output file

%% INPUT PARSING
pa=inputParser; pa.CaseSensitive = false; % Initialize input parser
addRequired(pa,'inputfile', @(x) ischar(x));
addRequired(pa,'xcalib', @(x) isreal(x) && isscalar(x));
addRequired(pa,'ycalib', @(x) isreal(x) && isscalar(x));
addParameter(pa,'outputfile',[], @(x) ischar(x));
parse(pa,inputfile,xcalib,ycalib,varargin{:}); % Read argument list
inputfile=pa.Results.inputfile;
xcalib=pa.Results.xcalib;
ycalib=pa.Results.ycalib;
outputfile=pa.Results.outputfile;

%% FIXED PARAMETERS OF INTERACTIVE PLOT
mks=20; % Marker size
smks=8; % Smaller marker size

%% INTERACTIVE DATA ACQUISITION IN WINDOW UNITS
%  in computer science coordinates such that:
%  [-] points are described by row-vectors;
%  [-] origin = upper left corner of window;
%  [-] ordinates are increasing downwards.

figure('NumberTitle','off','name',inputfile)
img=imread(inputfile); 
imshow(img); hold on;

button=0; i=0;
while (button~=3)
  [xx, yy, button] = ginput(1); % Capture only 1 mouse input from any button
  switch button
    case 1  % Add data point
      i=i+1; x(i)=xx; y(i)=yy;
      if (i<4)                             % Calibration points
                                           % inside blue squares
          plot(x(i),y(i),'sb',...          % First part of marker: square
               'markersize',mks); hold on;
          plot(x(i),y(i),'ob',...          % Second part of marker: 
               'markersize',smks,...       % filled circle inside square
               'markerfacecolor','g'); hold on;  
          
        
      else                                 % Data points inside red squares
	plot(x(i),y(i),'sr','markersize',mks); hold on; % Square
	plot(x(i),y(i),'or','markersize',smks,...       % Filled circle
             'markerfacecolor','g'); hold on;
      end
    case 2  % Remove previous data point
      if (i>0)
	marks=findobj(gca, 'type', 'line'); % Handles of markers in reverse capture order
	                                    % (lower index = last captured point)
	delete(marks(1:2));    % Two last objects deleted because marker has 2 compounds
	x(i)=[];y(i)=[];i=i-1; 
      end
    case 3           
      close; % Close imshow 
  end
end

%% CONVERSION TO ALGEBRA COORDINATES IN WINDOW UNITS
%  such that:
%  [-] points are described by column-vectors;
%  [-] origin = lower left corner of window;
%  [-] ordinates are increasing upwards.

R=[x' y']';     % Transforming to column vectors

R(2,:)=-R(2,:); % Mirroring to get the right ordinate axis orientation

origin = R(:,1);% The origin of the graph is the first digitized point
R=R-origin;     % Shift all data

%% CONVERSION TO ALGEBRA COORDINATES IN GRAPH UNITS

R(1,:)=R(1,:)*xcalib/R(1,2); % Abscissa calibration using abscissa of point 2
R(2,:)=R(2,:)*ycalib/R(2,3); % Ordinate calibration using ordinate of point 3

R(:,1:3)=[]; % Discard calibration points 1, 2 and 3

% Data points may be mouse-added randomly 
% because next command sorts the data
R=sortrows(R',1,'ascend')';              

%% SAVE RESULTS IN FILE
if isempty(outputfile)
    output_file=changefiletype(inputfile,'dat');
else
    output_file=outputfile;
end
printtable(R','file',output_file)% Row-wise save to file
fprintf('%d digitized data saved in file %s\n',size(R,2),output_file)

%% PLOT ACQUIRED DATA IN ORIGINAL GRAPH UNITS

figure('NumberTitle','off','name','Digitized curve')
plot(R(1,:),R(2,:)); hold on;
plot(R(1,:),R(2,:),'sr','markersize',smks); hold on;
axis auto
movegui(gcf,'north');

% End of function digitizer


 
