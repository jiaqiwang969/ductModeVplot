function sigma=variance(x)
% to caculate the variance of the signal
N=length(x);
m=meanvalue(x);
x=x-m;
sigma=sum(x.^2)/N;
sigma=sqrt(sigma);