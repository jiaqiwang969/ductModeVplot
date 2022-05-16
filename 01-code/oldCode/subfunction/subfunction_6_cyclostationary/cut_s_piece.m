function Suse=cut_s_piece(S,cutplace,fs,ShowStar,ShowEnd)

if nargin<4
    ShowStar=0;
    ShowEnd=floor(fs/2.56);
end


[r,c]=size(S);

cutpnum=round(cutplace*2*c/fs)+1;
Suse=S(:,cutpnum);

lend=length(Suse);
ShowStarN=floor(ShowStar*2*lend/fs)+1;
ShowEndN=floor(ShowEnd*2*lend/fs)+1;

dispSuse=Suse(ShowStarN:ShowEndN);
xld=(ShowStarN:ShowEndN+1)-1; 
xld=xld*fs/(2*lend);


if nargout==0
    titstr=['cyclic frequency = ',num2str(cutplace),' Hz'];
    figure(gcf+1);
    plot(xld,[abs(dispSuse') 0],'k');
    title(titstr);
    xlabel('frequency[Hz]');
    ylabel('amplitude[V]');
        axis tight;
    set(gcf,'color','white');
    
end

