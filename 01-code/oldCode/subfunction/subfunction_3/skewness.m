function sk=skewness(x)
% to caculate the skewness of the signal
N=length(x);
sigma=variance(x);
m=meanvalue(x);
x=(x-m)./sigma;
sk=sum(x.^3)/N;