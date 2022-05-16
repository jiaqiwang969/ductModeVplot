function [y]=create_gear_signal(N,fs,fz,fr,sign)

if nargin<5
    sign=0;
end

t=(0:N-1);
t=t/fs;
M=3;
L=5;
m=0:M;
l=0:L;
X=ones(1,M+1);
inph=0.1*pi*rand(1,M+1);
A=ones(1,L+1);
Aph=0.1*pi*rand(1,L+1);
B=ones(1,L+1);
Bph=0.1*pi*rand(1,L+1);

a=zeros(1,N);
b=zeros(1,N);

X=exp(-1*((0:M)*(-0.4*log2(0.1))/M));

A=5*exp(-1*((0:L)*(-0.2*log2(0.1))/L));

B=A;

if sign==0
    a=A*cos(2*pi*fr*l'*t+Aph'*ones(1,N));
elseif sign==1
    b=B*cos(2*pi*fr*l'*t+Bph'*ones(1,N));
else
    a=A*cos(2*pi*fr*l'*t+Aph'*ones(1,N));
    b=B*cos(2*pi*fr*l'*t+Bph'*ones(1,N));
end

y=(1+a).*(X*cos((2*pi*fz*m'*t+inph'*ones(1,N)+ones(M+1,1)*b)));


if nargout==0 
    figure(gcf+1);
    plot(t,y,'k');
    xlabel='time [S]';
    ylabel='amplitude [V]';
    title='the signal';
    showfft(y,fs);
    set (gcf,'Color','white');
end
