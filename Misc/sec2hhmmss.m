function [hh, mm, ss] = sec2hhmmss(t)

%% CONVERSION OF AMOUNT OF SECONDS TO HOURS:MINUTES:SECONDS
%
%% INPUT
%  t [second]
    
%% OUTPUT
%  hh [hour]
%  mm [minute]
%  ss [second]

if (~isreal(t))
    error('Wrong argument type.\n')
end
    
hh = floor(t / 3600);
t = t - hh * 3600;
mm = floor(t / 60);
ss = round(t - mm * 60);

end % End of function sec2hhmmss