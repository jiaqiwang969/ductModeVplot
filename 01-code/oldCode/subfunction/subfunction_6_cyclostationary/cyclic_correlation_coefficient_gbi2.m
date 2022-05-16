function [rho]=cyclic_correlation_coefficient_gbi2(x,Nt,Nhalftao,alpha,fs)

x_S0=x(1:Nt+Nhalftao*2);
x_S=x(1:Nt+Nhalftao*2);

y=showfft(x_S0,fs);
S0=abs(y).*abs(y);
S0_len=length(S0);
S0((S0<max(S0)/100))=0;

S=piece_analysis(x_S,alpha,fs,Nt,Nhalftao);
S=abs(S);
S_len=length(S);
S((S<max(S)/100))=0;

rho=zeros(1,S_len);
rhotmp=zeros(1,S_len);

for f=1:S_len
    fplus=(f-1)+round(alpha*2*S0_len/fs);
    if fplus>S0_len-1
        fplus=2*(S0_lenhalf-1)-fplus;
    end
    if S0(fplus+1)*S0(f)~=0 
        rho(f)=sqrt(S0(f)*S0(fplus+1));
    else
        rho(f)=0;
    end
end

if nargout==0
    len_rho=length(rho);
    figure(gcf+1);
    h=(0:len_rho-1)*fs/2/len_rho;
    plot(h,rho,'ko-');
    figure(gcf+1);plot(h,S);

end