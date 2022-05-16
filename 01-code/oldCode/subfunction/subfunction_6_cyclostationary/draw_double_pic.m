function draw_double_pic(x,fs,ShowStar,ShowEnd)

if nargin<3
    ShowStar=0;
    ShowEnd=fs/2.56;
end


lenx=length(x);
xd=x-mean(x);
x_star=0;           
x_end=0.5;
show_x_star=floor(x_star*fs)+1;    
show_x_end=floor(x_end*fs)+1;    
xshow=xd(show_x_star:show_x_end);
xshow_lnx=(1:length(xshow))-1;
xshow_lnx=xshow_lnx/fs;
xFFT=x(1:2048);
xtmp=xFFT-mean(xFFT);
xl=length(xFFT);
tmpy=fft(xtmp)/xl;
y=tmpy(1:floor(xl/2));
showy=y(1:floor(xl/2.56)+1);
xly=(1:floor(xl/2.56)+1)-1; 
xly=xly*fs/xl;
Nt=1024;
Nhalftao=512;
R=cyclic_cross_correlation_fast_gbi(x,x,Nt,Nhalftao);
d=sum(abs(R).^2)/sum(abs(R(:,1)).^2);
lend=length(d);
ShowStarN=floor(ShowStar*2*lend/fs)+1;
ShowEndN=floor(ShowEnd*2*lend/fs)+1;
dispd=d(ShowStarN:ShowEndN);
xld=(ShowStarN:ShowEndN+1)-1; 
xld=xld*fs/(2*lend);
if ShowStar==0
    dispd(1)=0;
    dispd(1)=max(dispd)+0.1*max(dispd);
end
S=scd(R);
[r,c]=size(S);  
dispS=S(1:r,1:c); 
xlS=(1:c+1)-1;              
xlS=xlS*fs/(2*c);        
ylS=(1:r+1)-1;  
ylS=ylS*fs/(2*r);
titletxt='Spectral Correlation Density';
xlabeltxt='cyclic frequency[Hz]';
ylabeltxt='frequency[Hz]';
zlabeltxt='CSD'; 
dispS=[dispS' zeros(c,1)];
dispS=[dispS' zeros(r+1,1)];      

figure(gcf+1);
set(gcf,'color','white'); 
subplot(2,1,1);
    plot(xshow_lnx,xshow,'k');
    title('the signal');
    xlabel('time[s]');
    ylabel('amplitude[V]');
    axis tight;
subplot(2,1,2);
    plot(xly,abs(showy),'k');
    title('The FFT of the Signal');
    xlabel('frequency[Hz]');
    ylabel('amplitude');
    axis tight;
figure(gcf+1);
set(gcf,'color','white'); 
subplot(2,1,1);
    contour(xlS,ylS,abs(dispS),90,'k');
    axis tight;
    set(gcf,'color','white');
    title(titletxt);
    xlabel(xlabeltxt);
    ylabel(ylabeltxt);
subplot(2,1,2);
    plot(xld,[dispd 0],'k');
    xlabel('cyclic frequency[Hz]');
    ylabel('DCS');
    title('degree of cyclostationarity');
    axis tight;

