function [feature]=feature_WE3(x)
[m,n]=size(x);
for i=1:n
    nwave=3;
    wpt1=wpdec(x(:,i),nwave,'db3'); %对数据进行小波包分解
    for j=1:2^nwave
        E(j)=norm(wpcoef(wpt1,[nwave,j-1]),2);%wpcoef(wpt1,[n,i-1])是求第n层第i个节点的系数
    end
    wave_p=E./sum(E);
    feature1=-wave_p.*log(wave_p);%        f_WE{i,k}=wave_p;
    feature(:,i)=[wave_p,feature1];
end
% feature=wave_p;
