function [m,i]=synchronous_average(x,fs,fsyn,numT)

fsynNtmp=round(fs/fsyn);
fsynN=fsynNtmp*numT;
M=floor(length(x)/fsynN);
xtmp=x(1:fsynN*M);
mtmp=zeros(M,fsynN);
mtmp=reshape(xtmp,[fsynN M]);
m=mean(mtmp');

i=M;
lnm=(1:fsynN)-1;
lnm=lnm*fs/length(x);
if nargout==0
    figure(gcf+1);
    plot(lnm,m);
    title('synchronous average');
    xlabel('time[s]');
    ylabel('amplitude');
end
