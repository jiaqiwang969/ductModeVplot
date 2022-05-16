function Rpiece=cyc_crocorr_certainalpha_gbi(x,alpha,fs,Nt,Ndt,Ntau)
%By Zhou (not sure)
%Ntau: length of timedelay
%Ndt


if nargin~=6
  error('Incorrect number of arguments for function cyclic_cross_correlation_fast');
end

[xrow,xcol]=size(x);
if ((xcol==0)|(xcol>2))
  error('x should be one or two column vector');
end

if fs/5.12<alpha 
    error('alpha is out of analysis frequency band.');
end

t=1:Ndt:Nt;
xrow=floor(xrow/Nt)*Nt;
tt=reshape(1:xrow,Nt,[]);
tt=tt(t,:);

lnt=length(t);
avg_num=floor(xrow/Nt);
Nlnt=lnt*avg_num;
R=zeros(Ntau,Nlnt);

for Rcol=1:Nlnt
    ti=tt(Rcol);
    taumax=min([ti-1,xrow-ti,ceil(Ntau/2)-1]);
    tau=-taumax:taumax;
    indices=rem(Ntau+tau,Ntau)+1;
    R(indices,Rcol)=x(ti+tau,1).*conj(x(ti-tau,xcol));   
    tau=ceil(Ntau/2);
    if ((ti<=xrow-tau)&(ti>=tau+1))
        R(tau+1,Rcol)=0.5*(x(ti+tau,1) * conj(x(ti-tau,xcol))  + x(ti-tau,1) * conj(x(ti+tau,xcol))) ; 
    end
end

RR=zeros(Ntau,lnt);
for M=1:avg_num
    tmpR=R(:,(M-1)*lnt+1:M*lnt);
    RR=RR+tmpR;
end
 RR=RR/avg_num;

meanR=mean(RR')';
meanR=meanR*ones(1,lnt);
RR=RR-meanR;

Rpiece=zeros(Ntau,1);
Palpha=round(alpha*Ndt*lnt/fs);                                   
Rpiece=exp(-i*pi*2*Palpha*(0:lnt-1)/lnt)*RR'/lnt;%

if nargout==0
    tau=-(ceil(Ntau/2)-1):floor(Ntau/2);
    fs_tau=fs/2;
    tau=tau/fs_tau;
    figure;
    plot(tau,abs(Rpiece),'k');
    set(gcf,'color','white');
    axis tight;
    figure;
    plot(tau,real(Rpiece),'k');
    set(gcf,'color','white');
    axis tight;
    figure;
    plot(tau,imag(Rpiece),'k');
    set(gcf,'color','white');
    axis tight;
end
    