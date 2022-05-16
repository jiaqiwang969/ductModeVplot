function [Spiece,pieceR]=cut_r_piece(R,cutplace,fs)


[r,c]=size(R);

cutpnum=round(cutplace*2*c/fs)+1;

pieceR=R(:,cutpnum);

Suse=fft(pieceR)/r;                        

Spiece=2*abs(Suse(1:floor(r/2)));

if nargout==0
    xl=(0:floor(r/2)-1)*fs/r;
    titstr=['alpha=',num2str(cutplace),'Hz'];
    figure('color','white');
    plot(xl,Spiece,'k');
  
    title(titstr);
    xlabel('frequency[Hz]');
    ylabel('amplitude');
    
    figure(gcf+1);
    plot((0:r-1)/fs,real(pieceR));
    title(titstr);
    xlabel('tao s');
    figure(gcf+1);
    plot((0:r-1)/fs,imag(pieceR));
    title(titstr);
    xlabel('tao s');
        figure(gcf+1);
    plot((0:r-1)/fs,abs(pieceR));
    title(titstr);
    xlabel('tao s');
end

