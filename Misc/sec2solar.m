function [dd,hh, mm, ss] = sec2solar(t)
%% CONVERSION OF AMOUNT OF SECONDS TO SOLAR TIME
%
%% INPUT
%  t [second]
    
%% OUTPUT
%  dd [day]
%  hh [hour]
%  mm [minute]
%  ss [second]
    
if (~isreal(t))
    error('Wrong argument type.')
end

hh = floor(t / 3600);
t = t - hh * 3600;
mm = floor(t / 60);
ss = round(t - mm * 60);

dd = floor(hh / 24);
hh = hh - dd*24;

end % End of function sec2solar