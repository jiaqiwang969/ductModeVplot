function y=show_cepstrum(x,fs,ShowStar,ShowEnd)

if nargin<3
    ShowStar=0;
    ShowEnd=fs/2.56;
end


xl=length(x);
xtmp=x-mean(x);
t1=fft(xtmp)/xl;
abst1=abs(t1);
abst1(abst1<10^-6)=10^-6;    
t2=20*log10(abst1);
ytmp=fft(t2-mean(t2))/length(t2);      
y=ytmp(1:length(ytmp));

ShowStarN=floor(ShowStar*xl/fs)+1;
ShowEndN=floor(ShowEnd*xl/fs)+1;
showy=y(ShowStarN:ShowEndN);
xly=(ShowStarN:ShowEndN)-1; 
xly=xly/fs;


if nargout==0
    figure(gcf+1);
    plot(xly,showy,'k');
    title('The Cepstrum of the Signal');
    xlabel('time[s]');
    ylabel('amplitude[dB]');
      axis tight;
    set(gcf,'color','white');  
end