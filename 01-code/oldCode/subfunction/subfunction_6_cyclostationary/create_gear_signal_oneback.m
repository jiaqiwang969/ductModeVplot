function [y,a,b]=create_gear_signal_one(N,fs,fz,fr,sign)

if nargin<5
    sign=0;
end

fra=fr(1);
frb=fr(length(fr));

% fra=fr;
% frb=fr;


t=(0:N-1);
t=t/fs;
L=5;
X=1;
l=0:L;
inph=0.1*pi;

A=ones(1,L+1);
Aph=0.3*pi*randn(1,L+1);
B=ones(1,L+1);
randn('seed',0);
Bph=pi/2+0.3*pi*randn(1,L+1);  

a=zeros(1,N);
b=zeros(1,N);

k=0:(80/L):80;

A=3*abs(randn(1))*exp(-0.03*k);
rand('seed',0);

randn('seed',0) ;
B=3*abs(randn(1))*exp(-0.03*k); 
rand('seed',0);



if sign==0
    a=A*cos(2*pi*fra*l'*t+Aph'*ones(1,N));
elseif sign==1
    b=B*cos(2*pi*frb*l'*t+Bph'*ones(1,N));
else
    a=A*cos(2*pi*fra*l'*t+Aph'*ones(1,N));
    b=B*cos(2*pi*frb*l'*t+Bph'*ones(1,N));
end


y=X*(1+a).*cos(2*pi*fz*t+inph+b);



if nargout==0 
    showfft(y,fs);
    showfft(a,fs);
    showfft(b,fs); 
end
