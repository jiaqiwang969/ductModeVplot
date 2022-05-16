function dispSpieces=draw_pieces(Spieces,fs,ShowStar,ShowEnd)

if nargin<3
    ShowStar=0;
    ShowEnd=fs/2.56;
end

[row col]=size(Spieces);

ShowStarN=floor(ShowStar*2*col/fs)+1;
ShowEndN=floor(ShowEnd*2*col/fs)+2;

dispSpieces=Spieces(:,ShowStarN:ShowEndN);
dispSpieces=20*log10(abs(dispSpieces));
xls=(ShowStarN:ShowEndN)-1; 
xls=xls*fs/(2*col);

txtstr=[];
for i=1:row;
    txtstr=[txtstr;'\leftarrow' num2str(i)];
end
figure;
plot(xls,dispSpieces);
text(zeros(1,row),dispSpieces(:,1),txtstr);
