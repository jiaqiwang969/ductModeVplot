function [y,a,b]=create_gear_signal_dgm(N,fs,fz,fr,sign)

if nargin<5
    sign=0;
end

t=(0:N-1);
t=t/fs;

% dt=randn(1,N-1)/4/fs; 
% dt=[0 dt];
% t=t+dt;

L=5;
X=1;
l=0:L;
inph=0.1*pi;

k=0:(20/L):20;
A=exp(-0.03*k);
B=exp(-0.03*k);

Aph=0.3*pi*randn(1,L+1);
Bph=pi/2+0.3*pi*randn(1,L+1);    

a=zeros(1,N);
b=zeros(1,N);



if sign==0
    a=A*cos(2*pi*fr*l'*t+Aph'*ones(1,N));
elseif sign==1
    b=B*cos(2*pi*fr*l'*t+Bph'*ones(1,N));
else
    a=A*cos(2*pi*fr*l'*t+Aph'*ones(1,N));
    b=B*cos(2*pi*fr*l'*t+Bph'*ones(1,N));
end

y=X*(1+a).*cos(2*pi*fz*t+inph+b);


if nargout==0 
    figure(gcf+1);
    plot(t,y,'k');
    xlabel='time [S]';
    ylabel='amplitude [V]';
    title='the signal';
    showfft(y,fs);
    set (gcf,'Color','white');
end
