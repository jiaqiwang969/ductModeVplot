function Rpiece=cyc_crocorr_certainalpha_gbi(x,alpha,fs,Nt,Ndt,Ntau)


if nargin~=6
  error('Incorrect number of arguments for function cyclic_cross_correlation_fast');
end

[xrow,xcol]=size(x);
if ((xcol==0)|(xcol>2))
  error('x should be one or two column vector');
end

t=1:Ndt:Nt;
xrow=2^(nextpow2(xrow)-1);
tt=1:Ndt:xrow;

lnt=length(t);
avg_num=floor(xrow/Nt);
Nlnt=lnt*avg_num;
R=zeros(Ntau,Nlnt);

for Rcol=1:Nlnt
    ti=tt(Rcol);
    taumax=min([ti-1,xrow-ti,round(Ntau/2)-1]);
    tau=-taumax:taumax;
    indices=rem(Ntau+tau,Ntau)+1;
    R(indices,Rcol)=x(ti+tau,1).*conj(x(ti-tau,xcol));
    tau=round(Ntau/2);
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

Rpiece=zeros(Ntau,1);
Palpha=floor(alpha*Nt/fs);
Rpiece=exp(-i*pi*2*Palpha*(0:lnt-1)/lnt)*RR'/lnt;

if nargout==0
    Rpiece=fftshift(Rpiece);
    tau=-(round(Ntau/2)-1):floor(Ntau/2);
    tau=tau/fs;
    figure(gcf+1);
    plot(tau,abs(Rpiece),'k');
    set(gcf,'color','white');
    axis tight;
end
    