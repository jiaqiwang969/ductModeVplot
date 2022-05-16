function Y=wvd_gbi(x,fs,Nt,halfNtao)

r=zeros(halfNtao*2+1,Nt);
for i=-halfNtao:halfNtao
    r(i+halfNtao+1,:)=x((i+halfNtao+1):(i+halfNtao+Nt)).*x(halfNtao+1:halfNtao+Nt);
end

R=fft(r)/(halfNtao*2+1);
Y=abs(R(1:halfNtao+1,:));

if nargout==0
    lt=0:Nt-1;
    lt=lt/fs;
    lf=0:halfNtao;
    lf=lf*fs/(halfNtao*2+1);
    figure(gcf+1);
    set(gcf,'color','white');
    contour(lt,lf,Y);
    axis tight;
    xlabel('time [s]');
    ylabel('frequency [Hz]');
    zlabel('amplitude [V]');
end
