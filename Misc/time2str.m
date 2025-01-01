function timestr = time2str(t)
%% CONVERSION OF AMOUNT OF SECONDS 
%% TO dd:hh:mm:ss STRING
    
%  Dependence on not built-in functions:
%  sec2solar.m

%% INPUT
%  t [seconds], real number

%% OUTPUT
%  timestr = 'dd:hh:mm:ss' [character string]
%  where
%   dd = number of days    : not printed if dd=0
%   hh = number of hours   : not printed if dd=hh=0
%   mm = number of minutes : not printed if dd=hh=mm=0
%   ss = number of seconds : always printed
    

if (~isreal(t) || ~isscalar(t))
    error('Argument of time2str(t) should be real scalar.')
end

if (t == Inf || t == NaN)
    timestr = sprintf('--:--:--:--');
else
    [dd,hh,mm,tt] = sec2solar(t);
    
    if dd==0 && hh==0 && mm==0 
        timestr=sprintf('%02d [s]',tt); return;
    end
    
    if dd==0 && hh==0 && mm>0
        timestr=sprintf('%02d:%02d [mm:ss]',mm,tt); return;
    end
    
    if dd==0 && hh>0 && mm>0
        timestr=sprintf('%02d:%02d:%02d [hh:mm:ss]',hh,mm,tt); return;
    end
    
    % if dd~=0, adjust width of 'days' field to the exact
    % required number of characters
    w=length(sprintf('%d',dd)); ds=[];
    if w==1; w=2; end % Esthetic symmetrization:
                      % Set width of 'days' field to 2
                      % if number of days < 10
    for i=1:w
        ds=[ds 'd'];
    end
    format=sprintf('%%0%dd:%%02d:%%02d:%%02d [%s:hh:mm:ss]',w,ds);
    timestr=sprintf(format,dd,hh,mm,tt);
end
end % End of function time2str
