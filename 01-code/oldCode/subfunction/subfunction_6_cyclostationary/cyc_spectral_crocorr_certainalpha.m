function [S,df]=cyc_spectral_crocorr_certainalpha(x,alpha,fs,g,M,L)



if ~exist('L')
    L=M;
end

if L<1 
    error('L must be no less than unity.');
end


[row col]=size(x);

if col>2
    x=x';
    [row col]=size(x);
end

lenx=row;
d_alpha=fs/lenx;
interval_f_N=round(alpha/d_alpha);
f_N=floor((lenx-interval_f_N-M)/L)+1;
lenx=(f_N-1)*L+M+interval_f_N; 
x=x(1:lenx,:);
X=fftshift(fft(x));
fnum=1:f_N;
m=(1:M)';
t=zeros(M,f_N);
t=m(:,ones(f_N,1))+(fnum(ones(M,1),:)-1)*L;


g=feval(g,M);
g=g/sum(g);
window_M=g(:,ones(f_N,1));
if col==1
    X1=zeros(M,f_N);X2=zeros(M,f_N);
    X1=X(t).*window_M;
    X2=X(t+interval_f_N).*window_M;
else
    X1=X(t,col);
    X1=reshape(X1,[M,f_N]);
    X1=X1.*window_M;
    X2=X(t+interval_f_N,1);
    X2=reshape(X2,[M,f_N]);
    X2=X2.*window_M;
end

if M==1
    S=conj(X1).*X2;
     S=S/lenx;
else
    Stmp=conj(X1).*X2;
    S=sum(Stmp);
     S=S/lenx;
end

if nargout==0
    df=-1*(fs/2-d_alpha*floor(M/2)-alpha/2):d_alpha*L:(fs/2-d_alpha*floor(M/2)-alpha/2);
    df=df(1:min(length(S),length(df)));
    S=S(1:min(length(S),length(df)));

    figure;
    subplot(2,1,1);
    plot(df,abs(S),'k');
    set(gcf,'color','white');
else
    df=-1*(fs/2-d_alpha*floor(M/2)-alpha/2):d_alpha*L:(fs/2-d_alpha*floor(M/2)-alpha/2);
    df=df(1:min(length(S),length(df)));
    S=S(1:min(length(S),length(df)));
end

    