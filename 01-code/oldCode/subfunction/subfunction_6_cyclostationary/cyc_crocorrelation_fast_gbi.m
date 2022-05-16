function [RR,alpha_decimation_factor]=cyc_crocorrelation_fast_gbi(x,Nt,Ndt,Ntau)


if nargin~=4
  error('Incorrect number of arguments for function cyclic_cross_correlation_fast');
end

[xrow,xcol]=size(x);
if ((xcol==0)|(xcol>2))
  error('x should be one or two column vector');
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

[Rrow,Rcol]=size(RR);

meanR=mean(R')';
meanR=meanR*ones(1,Rcol);

RR=fft(RR.').';
RR=RR/Rcol;


alpha_decimation_factor=1/Ndt;

if nargout==0  
    figure;
    fs=12000;
    fs_t=fs*alpha_decimation_factor;
    fs_tau=fs/2;
    tau=-(ceil(Ntau/2)-1):floor(Ntau/2);
    tau=fftshift(tau);
    surf((0:lnt-1)/lnt*fs_t,tau/fs_tau,abs(RR));
    set(gcf,'color','white');
end
