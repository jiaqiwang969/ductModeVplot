function drawRS(RS,d_ftau,d_alpha,sign,tysign)

if nargin<4
    sign='R';
    tysign='contour';
elseif nargin<5
    tysign='contour';
end

str_xlabel='cyclic frequency [Hz]';
if strcmp(sign,'R')
    str_title='cyclic autocorrelation';
    str_ylabel='tau [s]';
else
    str_title='spectral correlation density';
    str_ylabel='frequency [Hz]';
end

if strcmp(tysign,'plot3')
    [Rd_alpha,Pd_ftau]=meshgrid(d_alpha,d_ftau);
    str_draw='(Rd_alpha,Pd_ftau,RSshow)';
else
    str_draw='(d_alpha,d_ftau,RSshow,20)';
end

[RSrow,RScol]=size(RS);
Nftau=floor(RSrow/2)+1;
Na=floor(RScol/2)+1;
RSshow=abs(RS(1:Nftau,1:Na));   

figure;
subplot(2,1,1);
eval(strcat(tysign,str_draw));
set(gcf,'color','white');
title(str_title);
xlabel(str_xlabel);
ylabel(str_ylabel);
if strcmp(tysign,'contour')==0
    zlabel('amplitude [V]');
end

