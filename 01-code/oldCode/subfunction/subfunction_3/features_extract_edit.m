function Feature=features_edit(X)
%此次为paper2主要想实现功能为：包络能量
%首先需要均一化
%然后可能需要滤波
%
%par
%滤波器参数设置
% fs=204800
% f=[4*10/fs 6*20/fs 2*1000/fs 2*1020/fs];
% A=[0 1 0];
% rp=0.153;
% rs=16.92;
% devp=1-10^(-rp/20);
% devs=10^(-rs/20);
% dev=[devp devs devp];
% [n,f0,A0,w]=remezord(f,A,dev);
% if rem(n,2)
%    n=n+1;
% end
% b=remez(n,f0,A0,w);
% freqz(b,1,length(b),1);
% data_filter = filter(b,1,data);
load filter_fs204800_f40_1020.mat
x_filter = filter(b,1,X);
N = length(x_filter);
x = (x_filter - mean(x_filter))./std(x_filter);
Feature=Envelop2(x,204800);
end