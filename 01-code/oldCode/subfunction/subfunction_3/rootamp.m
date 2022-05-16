function xr=rootamp(x)
%to caculate the root amplitude of the signal
N=length(x);
xr=sum(sqrt(abs(x)))/N;
xr=xr.^2;