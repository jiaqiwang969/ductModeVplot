function y=cspectrum_gbi(x,fs)

x=x-mean(x);
N=length(x);
tmpy=abs(fft(x)/N);
tmpy(tmpy<10e-6)=10e-6;
tmpy2=ifft(20*log10(tmpy));
y=tmpy2(1:floor(length(tmpy2)/2)+1);

ly=1:length(y);
ly=(ly-1)/fs;
if nargout==0
    figure(gcf+1);
    plot(ly,abs(y));
end