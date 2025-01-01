function loop_progress(j,jmin,jmax,varargin)

%% VERY SIMPLE MONITORING OF THE PROGRESS OF A "FOR" LOOP
%% IN A TERMINAL WINDOW AS A ONE-LINE COMMAND
%  Updating on the {\em same} line: percentage of completed task,
%  estimated time to full completion and the averaged loop period.
%  Instead of the estimated time to completion, the final message
%  shows the loop duration.

%  NB: Because of the {\em wanted one-line command feature}, the
%  averaging of the loop period excludes the first iteration. For the
%  same reason, the final message shows the duration from the {\em
%  second} to the last iteration (See below: "Explanations"). 
%  Getting the total loop duration requires at least two commands: one
%  before and one after the loop: tic before the loop and toc
%  after the loop.

%% DEPENDENCE ON NOT BUILT-IN FUNCTIONS
%  time2str.m, sec2solar.m isintegervalue.m

%% INPUTS
%  j    = index of task progression (integer value)
%  jmin = index of first iteration  (integer value)
%  jmax = index of last  iteration  (integer value)
    
%% INSTRUCTIONS FOR USE
    
%  By interval of 1 second after the second iteration, this function
%  prints a message indicating the percentage of the task completion,
%  the estimated time to 100% completion and the averaged loop period.
%  Exceptions:
%  The task progression is never printed after the first iteration but
%  is always printed after the second iteration.
       
%  For a consistent timing, loop_progress must be invoked as last
%  instruction before the "end" statement of the monitored "for" loop:

%  for j=jmin:jmax
%      ...
%      ... loop's job ... 
%      ...
%      loop_progress(j,jmin,jmax)       
%  end 
    

%% EXPLANATIONS
%  Only the averaged loop period can be warranted because, in general,
%  the iteration time may vary for many of reasons (ex: due to
%  consequences of if-then-else implemented in an iteration).
    
%  To manage monitoring as a {\em one-line} command, it is not
%  possible to evaluate the time elapsed after the first iteration
%  because the clock is triggered {\em inside} the loop after this
%  first iteration. Therefore, the loop period can only be averaged
%  over iterations 2 to j, thereby {\em never} including the period of
%  the first iteration. This is a rather relevant feature since the
%  period of the first iteration is often not representative of the
%  period of the following iterations. Indeed the first iteration is
%  often much slower (for example due to some memory allocation that
%  will not be repeated in the next iterations).
    
%  The averaged loop period is used to update the estimated time to
%  completion. Therefore, if the iteration time decreases or increases
%  significantly as a function of the loop index j, the estimated time
%  to completion can not be a relevant approximation.

%% Recognized options in varargin 
% (uppercases for readability are optional): 
%
% if varargin{k} = 'PercentageStep', 
%                  then varargin{k+1} = real value [percent] 
%                  determining the frequency of loop progress
%                  messages after the second iteration instead of
%                  message every second.
%                  Ex: A value = 5 prints a message
%                      every time the task progressed by 5% 
%                      after the second iteration.
 
    
%% INPUT PARSING
pa=inputParser; pa.CaseSensitive = false;  
addRequired(pa,'j',    @(x) isintegervalue(x));
addRequired(pa,'jmin', @(x) isintegervalue(x));
addRequired(pa,'jmax', @(x) isintegervalue(x));
addParameter(pa,'percentagestep',[],...
             @(x) isreal(x) && isscalar(x) && x>0 && x<=100);
parse(pa,j,jmin,jmax,varargin{:}); 
pstep=pa.Results.percentagestep;

%% CORE JOB

% Values persisting until next call 
% to this function by the same script
persistent t0 period nb next

% Re-scaling as a "for i=1:imax" loop
if (jmin<=jmax)
    imax=jmax-jmin+1; i=j-jmin+1; 
else
    imax=jmin-jmax+1; i=j-jmax+1;
end

% From here, using re-scaled indexes i and imax

per = 100*i/imax; % Current percentage of task completion

t = clock(); % Get current time

if (i==1)     
    t0 = t;    % Since loop_progress must be invoked before the "end"
               % statement of the monitored "for" loop, t0 is here
               % the time after the first iteration was completed.
    period = 0; dur= Inf;                  % Initialization
    msg=sprintf([' %05.2f%% | Est. time to completion: %s | Av. loop ' ...
                 'period = %#12.4G [s]'],per,time2str(dur),period);
    fprintf('%s',msg); % First printed message
    nb=length(msg);    % Needed to move cursor back to start of
                       % line before printing next message
    return;
end

if (i==2) % i=2 always prints first value of period
    period=etime(clock(),t0); % Time of iteration i=2
    dur = period * (imax-i);  % Est. time to completion  
    for j=1:nb
        fprintf('\b');        % Move cursor back to start of line
    end
    msg=sprintf([' %05.2f%% | Est. time to completion: %s | Av. loop ' ...
                 'period = %#12.4G [s]'],per,time2str(dur),period);
    fprintf('%s',msg);        % Update message on same line
    nb=length(msg);           % Needed to move cursor back to start of
                              % line before printing next message 
    if isempty(pstep)
        next=2+round(1/period); % Next percentage to be printed
    else
        next=ceil(pstep*ceil(per/pstep)*imax/100);
    end
    return;
end

if (i==next)  
    period=etime(clock(),t0)/(i-1); % Averaged loop period  
    dur = period * (imax-i);        % Est. time to completion
    for j=1:nb
        fprintf('\b');        % Move cursor back to start of line
    end
    msg=sprintf([' %05.2f%% | Est. time to completion: %s | Av. loop ' ...
                 'period = %#12.4G [s]'],per,time2str(dur),period);
    fprintf('%s',msg); % Update message on same line
    nb=length(msg);    % Needed to move cursor back to start of
                       % line before printing next message
    if isempty(pstep)
        next=next+round(1/period);% Next index i to be printed
    else
        next=ceil(pstep*ceil(per/pstep)*imax/100);
    end
    return;
end

if (i>=imax) % Message when task is completed
    for j=1:nb
        fprintf('\b') % Move cursor back to start of line
    end  
    fprintf([' %06.2f%% '...
             '| Loop duration   = %s (2d to last iteration)\n'],...
            100,time2str(etime(t,t0)));
    fprintf(['         '...
             '| Av. loop period = %#12.4G [s] '...
             '(excl. 1st iteration)\n\n'],...
            period);
        
    return;
end

end % End of function loop_progress