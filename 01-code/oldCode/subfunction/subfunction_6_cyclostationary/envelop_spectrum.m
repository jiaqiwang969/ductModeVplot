% envelop
function [X]=envelop_spectrum(x,fs,ShowStar,ShowEnd)

if nargin<3
    ShowStar=0;
    ShowEnd=floor(fs/2.56);
end

[rol col]=size(x);
if rol>col
    x=x';
end

x=x-mean(x);

Xh=abs(hilbert(x));
lnX=length(Xh);
X=fft(Xh-mean(Xh))/lnX;
X=abs(X);
Xshow=X;

ShowStarN=floor(ShowStar*lnX/fs)+1;
ShowEndN=floor(ShowEnd*lnX/fs)+1;
Xdisp=Xshow(ShowStarN:ShowEndN);
xld=(ShowStarN:ShowEndN)-1; 
xld=xld*fs/lnX;


if nargout==0
      figure;
   set(gcf,'Position',[364   456   672   252]);axes( 'position' , [0.11,0.2,0.82,0.65] , 'box' , 'on');
    plot(xld,Xdisp,'k');
    set(gcf,'color','white');
    set(gca,'fontsize',12);
    set(get(gca,'xlabel'),'string','\fontsize{12}f (Hz)');
    set(get(gca,'ylabel'),'string','\fontsize{12}Envelope Spectrum');
end
    