function Y=short_fft(x,fs,win_len,shift_len,win_type,f_len)


lx=length(x);
dt_len=floor((lx-win_len)/shift_len);
df_len=win_len;
dt=dt_len/fs;
df=fs/df_len;

X=zeros(df_len,dt_len);
for n=0:dt_len-1
    X(:,n+1)=x(1+n*shift_len:(n*shift_len+df_len))';
end

win=eval([win_type,'(',num2str(win_len),')']);
win=win(:,ones(1,dt_len));
X=X.*win;

Y=fft(X,f_len)/f_len;

if nargout==0
    surf((1:dt_len)*dt,(1:floor(f_len/2)+1)*fs/f_len,abs(Y));
    set(gcf,'color','white');
    set(gca,'fontsize',12);

   xlabel('\fontsize{12}time [s]');
    ylabel('\fontsize{12}frequency [Hz]');
    zlabel('\fontsize{12}amplitude [v]');
    axis tight;
 
end

 