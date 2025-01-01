function plotzeros(varargin) 

%% TRACING X=0 AND Y=0 AXES ON CURRENT PLOT
      
%  varargin may contain line properties (color, width, style)
    
if isempty(varargin)
    varargin={'color',[0 0 0],'linewidth',0.5,'linestyle','-'};
end
    
h = findall(0,'type','figure');
if isempty(h)
    error('function plotzeros: can not draw on an empty plot')
end

xline = @(xval, varargin) line([xval xval], ylim, varargin{:});
yline = @(yval, varargin) line(xlim, [yval yval], varargin{:});
xline(0,varargin{:});
yline(0,varargin{:});

end % End of function plotzeros