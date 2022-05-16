function ku=kurtosis(x)
% to caculate the kurtosis of the signal
N=length(x);
sigma=variance(x);
m=meanvalue(x);
x=(x-m)./sigma;
ku=sum(x.^4)/N;