mystartdefaults
func = @(x) exp(-0.2*x).*cos(x);
x= 0:0.01:10; n = length(x);
y = func(x);
noise = rand(1,n)-0.5;

A=1;
yref = y+A*noise;

stats = gof(y,yref, 'short', true)