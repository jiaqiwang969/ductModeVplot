function xrms=rootmeansquare(x)
%to caculate the root mean square of the signal
N=length(x);
xrms=sum(x.^2)/N;
xrms=sqrt(xrms);