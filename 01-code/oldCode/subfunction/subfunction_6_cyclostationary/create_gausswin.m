function g=create_gausswin(N,evn,avr)

sign=0;
if nargin<2
    evn=0;
    avr=1/(2*pi);
    sign=1;
end


half_w_len=floor(N/2)+1;            
nhalf_w_len=ceil(N/2)-1;            


n=0:half_w_len-1;
gt=n*(3*sqrt(avr)/half_w_len);
g=zeros(1,N);
tmpg=(1/sqrt(2*pi*avr))*exp(-0.5/avr*(gt-evn).*(gt-evn));

if sign==1
    tmpg=tmpg*sqrt(sqrt(2));                   
end

tmpg_ver=tmpg(half_w_len:-1:1);
g=[tmpg_ver(1:nhalf_w_len),tmpg];
len_g=length(g);


if nargout==0
    figure(gcf+1);
    set(gcf,'color','white');
    plot((-1*nhalf_w_len:(half_w_len-1)),g,'o-');
end
