function [feature]=feature_WE3(x)
[m,n]=size(x);
for i=1:n
    nwave=3;
    wpt1=wpdec(x(:,i),nwave,'db3'); %�����ݽ���С�����ֽ�
    for j=1:2^nwave
        E(j)=norm(wpcoef(wpt1,[nwave,j-1]),2);%wpcoef(wpt1,[n,i-1])�����n���i���ڵ��ϵ��
    end
    wave_p=E./sum(E);
    feature1=-wave_p.*log(wave_p);%        f_WE{i,k}=wave_p;
    feature(:,i)=[wave_p,feature1];
end
% feature=wave_p;
