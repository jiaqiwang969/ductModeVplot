function [rho,t]=cyclic_correlation_coefficient_gbi(x,Nt,Nhalftao,multiNhalftao,alpha,fs)
% can't work!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

S0_len=multiNhalftao*Nhalftao+1;
S0_len=512;
x_S0=x(1:S0_len);
y=showfft(x_S0,fs);
S0=abs(y).*abs(y);
S0_lenhalf=length(S0);
S0((S0<max(S0)/100))=0;
% figure(gcf+1);plot(S0);

x_S=x(1:Nt+Nhalftao*2);
S=piece_analysis(x_S,alpha,fs,Nt,Nhalftao);

S=abs(S);
S_len=length(S);

S((S<max(S)/100))=0;

t=zeros(1,S_len);
for f=1:S_len
    if abs((f-1)*fs/S_len/2-85)<0.1
        tmpf=(f-1)*fs/S_len/2;
        tmpf=tmpf;
    end
    if abs((f-1)*fs/S_len/2-95)<0.1
        tmpf=(f-1)*fs/S_len/2;
        tmpf=tmpf;
    end
    
    fplus= round((f-1)*S0_lenhalf/S_len)+round(alpha*S0_len/fs);
    fminus=round((f-1)*S0_lenhalf/S_len)-round(alpha*S0_len/fs);
    
    if fplus>S0_lenhalf-1
        fplus=2*(S0_lenhalf-1)-fplus;
    end
    if fminus<0
        fminus=-1*fminus;
    end
 
    if S0(fplus+1)*S0(fminus+1)~=0 
        rho(f)=S(f)/sqrt(S0(fplus+1)*S0(fminus+1));
        t(f)=sqrt(S0(fplus+1)*S0(fminus+1));
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