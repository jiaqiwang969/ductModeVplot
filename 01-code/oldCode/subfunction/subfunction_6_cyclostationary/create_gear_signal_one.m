function [y,a,b]=create_gear_signal_one(N,fs,fz,fr,sign)

if nargin<5
    sign=0;
end

fra=fr(1);frb=fr(length(fr));
C=1;

amp=[0.1 0.1];amp_a=amp(1);amp_b=amp(length(amp));      
L=[3 3];La=L(1);Lb=L(length(L));      
X=1;inph=0.1*pi;     


t=(0:N-1)/fs; la=0:La;lb=0:Lb;
A=ones(1,La+1);
Aph=0.3*pi*randn(1,La+1);
B=ones(1,Lb+1);
Bph=pi/3+0.3*pi*randn(1,Lb+1);   
a=zeros(1,N);b=zeros(1,N);

k=0:(80/La):80;
A=amp_a*exp(-0.03*k);
A=[A(length(A)) A(1:length(A)-1)];              
k=0:(80/Lb):80;
B=amp_b*exp(-0.03*k); 
B=[B(length(B)) B(1:length(B)-1)];

A=[0.1 0.7 0.4 0.3]*1; 
Aph=[0.5 0.3 0.05 0]*pi;
B=[0.1 0.8 0.4 0.2]*1;     
Bph=[0.2 0.1 0.4 0.08]*pi+pi/3;

if sign==0
    a=A*cos(2*pi*fra*la'*t+Aph'*ones(1,N));
elseif sign==1
    b=B*cos(2*pi*frb*lb'*t+Bph'*ones(1,N));
else
    a=A*cos(2*pi*fra*la'*t+Aph'*ones(1,N));
    b=B*cos(2*pi*frb*lb'*t+Bph'*ones(1,N));
end


y=X*(C+a).*cos(2*pi*fz*t+inph+b);



if nargout==0 
    figure;
    plot(t,y,'k');
    xlabel='time [S]';
    ylabel='amplitude [V]';
    title='the signal';
    showfft(y,fs);
end
