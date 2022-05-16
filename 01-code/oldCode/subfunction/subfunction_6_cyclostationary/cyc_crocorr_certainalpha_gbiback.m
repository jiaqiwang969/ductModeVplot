function Rpiece=cyc_crocorr_certainalpha_gbi(x,alpha,fs,Nt,Ntau)


if nargin~=5
  error('Incorrect number of arguments for function cyclic_cross_correlation_fast');
end

[xrow,xcol]=size(x);
if ((xcol==0)|(xcol>2))
  error('x should be one or two column vector');
end

NNt=floor(xrow/Nt)*Nt;
R=zeros(Ntau,NNt);

for Rcol=1:NNt
    taumax=min([Rcol-1,xrow-Rcol,round(Ntau/2)-1]);
    tau=-taumax:taumax;
    indices=rem(Ntau+tau,Ntau)+1;
    R(indices,Rcol)=x(Rcol+tau,1).*conj(x(Rcol-tau,xcol));
    tau=round(Ntau/2);
    if ((Rcol<=xrow-tau)&(Rcol>=tau+1))
        R(tau+1,Rcol)=0.5*(x(Rcol+tau,1) * conj(x(Rcol-tau,xcol))  + x(Rcol-tau,1) * conj(x(Rcol+tau,xcol))) ;
    end
end
RR=zeros(Ntau,Nt);
for M=1:floor(xrow/Nt)
    tmpR=R(:,(M-1)*Nt+1:M*Nt);
    RR=RR+tmpR;
end
RR=RR/floor(xrow/Nt);

Rpiece=zeros(Ntau,1);
Palpha=floor(alpha*Nt/fs);
Rpiece=exp(-i*pi*2*Palpha*(0:Nt-1)/Nt)*RR'/Nt;

if nargout==0
    Rpiece=fftshift(Rpiece);
    tau=-(round(Ntau/2)-1):floor(Ntau/2);
    tau=tau/fs;
    figure(gcf+1);
    plot(tau,abs(Rpiece),'k');
    set(gcf,'color','white');
    axis tight;
end
    