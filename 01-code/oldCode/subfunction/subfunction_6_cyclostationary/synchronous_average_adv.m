function [m,i]=synchronous_average_adv(x,fs,fsyn,numT)

x=x-mean(x);
[col row]=size(x);
if col>row
    x=x';
end

dt=1/fs;
T=1/fsyn*numT;
N=floor(T/dt);
lnx=length(x);

M=1;
i=1;
synx=[];
while lnx>M+N
    tmpN=floor((M-1)*dt/T);
    if (i*T-(M-1+N)*dt)/dt>0.5
        M=M+1;
    elseif (i*T-(M-1+N)*dt)/dt<-0.5
        M=M-1; 
    end
        tmpx=x(M:M+N-1);
    try
        synx=[synx; tmpx];
    catch
        lentmp=length(tmpx);
    end
        M=M+N;
        i=i+1;
end 
i=i;
if i==2
    m=synx;
else
    m=mean(synx);
end

lnm=(1:N)-1;
lnm=lnm*fs/length(x);
if nargout==0
    figure(gcf+1);
    set(gcf,'color','white');
    plot(lnm,m,'b');
    title('synchronous average');
    xlabel('time[s]');
    ylabel('amplitude');
end
