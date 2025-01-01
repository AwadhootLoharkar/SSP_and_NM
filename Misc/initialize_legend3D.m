function initialize_legend3D()
%%% Initialize legend in current 3D figure
%   by creating an invisible graphical object
%   inside the current figure window.
%   The 'DisplayName' of this invisible object
%   is then used to print the title (here 'Legend')
%   of the legend box as first line in this box.

    plot3(0,0,0,'o','color','w','DisplayName','\bf Legend'); 
    hold on; 
end
