mystartdefaults

%% fp.m : EXPLORING CRITICAL VALUES OF FLOATING-POINT NUMBERS

graph_fp_density=true;

%% INPUT DATA 

%% Single precision
q(1) = 8;   % Number of bits for exponent 
p(1) = 23;  % Number of bits for significand

%% Double precision 
q(2) = 11;  % Number of bits for exponent 
p(2) = 52;  % Number of bits for significand

for i=1:2
    if(i==1)
        fprintf('\nSINGLE PRECISION\n')
    else
        fprintf('\nDOUBLE PRECISION\n')
    end
    fprintf('2^(%d+1)-1 possible values of mantissa = %d\n',p(i),2^(p(i)+1)-1);% Add phantom & use eq.(3.33)
    sd(i)=log10(2^(p(i)+1)-1); % Unfortunately, not an integer power of 10
    fprintf('log10 of above value = significant decimal digits = %d\n',sd(i)); 
    fprintf('In practice: %d significant decimal digits, %dth decimal digit uncertain\n',floor(sd(i)),floor(sd(i))+1);
    fprintf('Displaying more than %d significant decimal digits is meaningless\n',floor(sd(i))+1)
end

%% MACHINE EPSILON

fprintf('\nMACHINE EPSILON PREDEFINED IN MATLAB/OCTAVE\n');
fprintf('eq.(3.12) for single precision              = %15.7E\n',single(2)^(-single(p(1))));
fprintf('Intrinsic function eps for single precision = %15.7E\n',eps(single(1.0)));			  
fprintf('\neq.(3.12) for double precision              = %23.15E\n',2^(-p(2)));    
fprintf('Intrinsic function eps for double precision = %23.15E\n',eps(1.0));

%%  TABLE OF CRITICAL VALUES

for i=1:2  
    eps_0(i)=2^(-p(i));                % eq. (3.12) machine epsilon
    excess = 2^(q(i)-1)-1;             % Excess for offset binary - 1 representation (p. 41)
    nmax(i)= 2^(q(i))-1 - excess - 1 ; % eq.(3.25) Maximum non special value of the exponent
    nmin(i)= 1 - excess;               % eq.(3.25) Minimum non special value of the exponent
    v(i)= (2-eps_0(i)) * 2^nmax(i);    % eq.(3.28) Overflow level
    u(i)= 2^nmin(i);                   % eq.(3.29) Underflow level
    eps_nmin(i)= 2^(nmin(i)-p(i));     % eq.(3.30) Smallest positive subnormal number
end

header{1}='Precision'; header{2}='nmin'; header{3}='nmax';
header{4}='eps_0';     header{5}='v';    header{6}='u';   header{7}='eps_nmin';
printtable( [(1:2)' nmin' nmax' eps_0' v' u' eps_nmin'],...
            'precision',3,...
            'integer',[1 2 3],...
            'colname',header,...
            'title','TABLE OF CRITICAL VALUES')


fprintf('Double precision: u - realmin = %23.15E\n',u(2)-realmin);
fprintf('Double precision: v - realmax = %23.15E\n',v(2)-realmax);

%% SUBNORMAL NUMBERS

fprintf('\nDOUBLE PRECISION: SUBNORMAL NUMBERS\n');
fprintf(' k    u/(2^k)                 significant      bit pattern\n');
fprintf('                            decimal digits\n');
for k=0:63
    sd(k+1)=floor(log10(2^(p(2)-k))); % significant decimal digits
    if(sd(k+1)<0)
        sd(k+1)=0;
    end
    fprintf('%2d  %23.15E %11d      %s\n',k,u(2)/(2^k),sd(k+1),float2bin(u(2)/(2^k))); 
end

%% INTERVALS BETWEEN FLOATING-POINT NUMBERS & DENSITY

kmax = 308;
for k=-kmax:1:kmax
    i=kmax+k+1;
    x=10^k; xx(i)=x;
    if (x <= u(2))
        x = u(2);
    end
    vk(i)=k;
    epsilon(i,1) = 2^(floor(log2(x))-p(2)); % eq.(3.11) combined with eq.(3.10)
    epsilon(i,2) = eps(x);                  % Intrinsic function in Matlab/Octave = epsilon(:,1)
    epsilon(i,3) = eps_0(2)*abs(x);         % eq. (3.18) = majoration of epsilon(:,1)
    density(i)   =-log10(epsilon(i,1));     % Log10 of density of
                                            % floating-point numbers per unit length
end

header{1}='x';header{2}='epsilon(x)'; header{3}='eps(x)';
header{4}='abs(x)*eps_0';header{5}='log10(density)';
printtable([xx' epsilon density'],...
           'precision',4,...
           'normalized',true,...
           'colname',header,...
           'title','Double precision: Intervals between floating-point numbers & density')

if graph_fp_density
    plot(vk,density);
    xlim([-kmax,kmax]);
    xlabel('Power of 10');
    ylabel('log10 of floating-point number density');
    hold on;
end

%% Inf & NaN

fprintf('\nHOW Inf ARISES\n');
fprintf('2^(nmax+1)=2^%4d in single = %d ; bit pattern = %s\n',...
	single(nmax(1)+1),2^(single(nmax(1)+1)),float2bin(2^(single(nmax(1)+1))));
fprintf('2^(nmax+1)=2^%4d in double = %d ; bit pattern = %s\n',...
	nmax(2)+1,2^(nmax(2)+1),float2bin(2^(nmax(2)+1)));

fprintf('\nWhat if exponent of base=2 grows above nmax+1 ?\n');
fprintf('2^(nmax+2)=2^%4d in single = %d\n',single(nmax(1)+2),2^(single(nmax(1)+2)));
fprintf('2^(nmax+2)=2^%4d in double = %d\n',nmax(2)+2,2^(nmax(2)+2));

fprintf('\nEXAMPLES OF NaN EVENTS IN DOUBLE PRECISION\n');
fprintf('O/O      = %d ; bit pattern = %s\n',0/0,float2bin(0/0));
fprintf('0 * Inf  = %d\n',0*2^(nmax(2)+2));
fprintf('Inf / Inf= %d\n',2^(nmax(2)+1)/2^(nmax(2)+2));
fprintf('Inf - Inf= %d\n',2^(nmax(2)+1)-2^(nmax(2)+2));

%% EXAMPLE OF Inf EVENTS & LOSS OF ACCURACY WHEN USING EXPONENTIAL FUNCTION

fprintf('\nDOUBLE PRECISION: EXAMPLE OF Inf EVENTS & LOSS OF ACCURACY WHEN USING EXPONENTIAL FUNCTION\n');

argmax=log(v(2));
fprintf('\nLogarithm of double precision overflow level = %23.15E\n',argmax);
fprintf('epsilon(%23.15E)=%23.15E\n',argmax,eps(argmax));

j=0;
for i=0:5  
    j=j+1;
    a(j)=argmax;
    delta(j)=eps(argmax)/(2^i);
    x(j)=a(j)+delta(j);
    y(j)=exp(x(j));
end
for i=5:-1:-5 
    j=j+1;
    a(j)=argmax;
    delta(j)=-eps(argmax)/(2^i);
    x(j)=a(j)+delta(j);
    y(j)=exp(x(j));
end

header{1}='argmax'; header{2}='delta'; 
header{3}='x=argmax+delta'; header{4}='exp(x)'; 
printtable([a' delta' x' y'],...
           'precision',15,...
           'normalized',true,...
           'colname',header)

%% "OFFSET BINARY - 1" prefered in IEEE 754 standard
%   because 1/u does not overlow (1/u<v) and 1/v < u (subnormal number)
fprintf('\nOFFSET BINARY -1\n')

for i=1:2
    eps_0(i)=2^(-p(i));                % eq. (3.12) machine epsilon
    excess = 2^(q(i)-1)-1;             % Excess for offset binary - 1 representation (p. 33)
    nmax(i)= 2^(q(i))-1 - excess - 1 ; % eq.(3.25) Maximum non special value of the exponent
    nmin(i)= 1 - excess;               % eq.(3.25) Minimum non special value of the exponent
    v(i)= (2-eps_0(i)) * 2^nmax(i);    % eq.(3.28) Overflow level
    u(i)= 2^nmin(i);                   % eq.(3.29) Underflow level
    fprintf('\nPrecision   nmin   nmax \n');
    fprintf('%9d   %4d  %4d  | u    =%23.15E | 1/u  =%23.15E \n',i,nmin(i),nmax(i),u(i),1/u(i));
    fprintf('%9d   %4d  %4d  | 1/v  =%23.15E | v    =%23.15E \n',i,nmin(i),nmax(i),1/v(i),v(i));
    fprintf('%9d   %4d  %4d  | 1/u-v=%23.15E | 1/v-u=%23.15E \n',i,nmin(i),nmax(i),1/u(i)-v(i),1/v(i)-u(i) );
end

%% "OFFSET BINARY" discarded in IEEE 754 standard
%   because 1/u overflows (1/u>v) and 1/v <= u

fprintf('\nOFFSET BINARY\n')
for i=1:2
    eps_0(i)=2^(-p(i));                % eq. (3.12) machine epsilon
    excess = 2^(q(i)-1);               % Excess for offset binary representation (p. 33)
    nmax(i)= 2^(q(i))-1 - excess - 1 ; % eq.(3.25) Maximum non special value of the exponent
    nmin(i)= 1 - excess;               % eq.(3.25) Minimum non special value of the exponent
    v(i)= (2-eps_0(i)) * 2^nmax(i);    % eq.(3.28) Overflow level
    u(i)= 2^nmin(i);                   % eq.(3.29) Underflow level
    fprintf('\nPrecision   nmin   nmax \n');
    fprintf('%9d   %4d  %4d  | u    =%23.15E | 1/u  =%23.15E \n',i,nmin(i),nmax(i),u(i),1/u(i));
    fprintf('%9d   %4d  %4d  | 1/v  =%23.15E | v    =%23.15E \n',i,nmin(i),nmax(i),1/v(i),v(i));
    fprintf('%9d   %4d  %4d  | 1/u-v=%23.15E | 1/v-u=%23.15E \n',i,nmin(i),nmax(i),1/u(i)-v(i),1/v(i)-u(i) );
end

