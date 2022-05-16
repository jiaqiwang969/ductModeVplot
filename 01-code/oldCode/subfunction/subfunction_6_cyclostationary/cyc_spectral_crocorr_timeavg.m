function S=cyc_spectral_crocorr_timeavg(x,fs,f_N,a,g,L)


if ~exist('L')
  L=1;
end
if floor(f_N/2)*2~=f_N
    warning('f_N must be odd.');
    f_N=f_N-1;
end
[row col]=size(x);

if col>2
    x=x';
    [row col]=size(x);
end

lenx=row;
alpha_N=floor((lenx-f_N)/L)+1;
lenx=(alpha_N-1)*L+f_N;
x=x(1:lenx,:);
a=feval(a,f_N);
g=feval(g,alpha_N);
a=a/sum(a);
win_a=a(:,ones(alpha_N,1));
g=g'/sum(g);


alphanum=1:alpha_N;
fn=(1:f_N)';
t=zeros(f_N,alpha_N);
t=fn(:,ones(alpha_N,1))+(alphanum(ones(f_N,1),:)-1)*L;

if col==1
    x1=zeros(f_N,alpha_N);x2=zeros(f_N,alpha_N);
    x1=x(t).*win_a;
    x2=x1;
else
    x1=x(t,1);
    x1=reshape(x1,[f_N,alpha_N]);
    x1=x1.*win_a;
    x2=x(t,col);
    x2=reshape(x2,[f_N,alpha_N]);
    x2=x2.*win_a;
end

X1=fftshift(fft(x1));
X2=fftshift(fft(x2));


S=zeros(2*f_N-1,(2*f_N-1));

for k=(-f_N+1):0
    i=k+f_N;
    tmpX1=X1(1:i,:);
    tmpX2=X2((f_N-i+1):f_N,:);
    win_g=g(ones(i,1),:);
    tmpS=tmpX1.*conj(tmpX2).*win_g;
    S(f_N-i+1:2:f_N+i-1,i)=sum(tmpS')';
    clear tmpS;
    clear tmpX1;
    clear tmpX2;
    clear win_g;
end

for k=1:(f_N-1)
    i=k+f_N;
    tmpX1=X1((i-f_N+1):f_N,:);
    tmpX2=X2(1:f_N-i+f_N,:);
    win_g=g(ones(f_N-i+f_N,1),:);
    tmpS=tmpX1.*conj(tmpX2).*win_g;
    tmp=i-f_N+1:2:i+f_N-2;
    S(i-f_N+1:2:2*f_N-(i-f_N)-1,i)=sum(tmpS')';
    clear tmpS;
    clear tmpX1;
    clear tmpX2;
    clear win_g;
end


if nargout==0
    xa=-1*(f_N-1):(f_N-1);
    xa=xa*fs/length(xa);
    ya=-1*(f_N-1):(f_N-1);
    ya=ya*2*fs/length(ya);
    [Yl,Xl]=meshgrid(ya,xa);  
    figure;
    contour(Xl,Yl,abs(S));
end
