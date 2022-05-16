function y=showfft(x,fs,ShowStar,ShowEnd)

if nargin==0
    [fname, pname] = uigetfile( 'E:\GBi\matlab\data\geardata\helicalgear4137_mat\*.mat', 'Open a file');
    fp=[pname, fname];
    load(fp);
    clear fname pname fp;
    fs=1;
    ShowStar=0;
    ShowEnd=fs/2.56;
elseif nargin<2
    fs=1;
    ShowStar=0;
    ShowEnd=fs/2.56;
elseif nargin<3
    ShowStar=0;
    ShowEnd=fs/2.56;
end


xl=length(x);
xtmp=x;
xtmp=x-mean(x);
tmpy=fft(xtmp)/xl;
y=tmpy;



if nargout==0
    ShowStarN=floor(ShowStar*xl/fs)+1;
    ShowEndN=floor(ShowEnd*xl/fs)+1;
    showy=y(ShowStarN:ShowEndN);
    xly=(ShowStarN:ShowEndN)-1; 
    xly=xly*fs/xl;

   figure;
   subplot(2,1,1);
    plot(xly,abs(showy),'k');
    title('\fontsize{12}the maplitude spcetrum');
    xlabel('\fontsize{12}frequency[Hz]');
    ylabel('\fontsize{12}amplitude[v]');
      axis tight;
    set(gcf,'color','white');  
    set(gca,'fontsize',12);
end