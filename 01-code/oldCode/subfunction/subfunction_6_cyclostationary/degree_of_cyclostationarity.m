function d=degree_of_cyclostationarity(R,fs,ShowStar,ShowEnd)

if nargin<3
    ShowStar=0;
    ShowEnd=floor(fs/2.56);
end


d=sum(abs(R).^2);


if nargout==0
    
lend=length(d);
ShowStarN=floor(ShowStar*2*lend/fs)+1;
ShowEndN=floor(ShowEnd*2*(lend-1)/fs)+1;   

dispd=d(ShowStarN:ShowEndN);
xld=(ShowStarN:ShowEndN)-1; 
xld=xld*fs/(2*(lend-1));

d=d(ShowStarN:ShowEndN);

if ShowStar==0
    dispd(1)=0;
    dispd(1)=max(dispd)+0.1*max(dispd);
end


    figure(gcf+1);
    plot(xld,dispd,'k');
    xlabel('\fontsize{12}cyclic frequency[Hz]');
    ylabel('DCS');
    title('\fontsize{12}degree of cyclostationarity');
    axis tight;
    set(gcf,'color','white');
        set(gca,'fontsize',12);

end