function mamp=meanamp(x)
%to caculate the mean amplitude of the signal
N=length(x);
mamp=sum(abs(x))/N;