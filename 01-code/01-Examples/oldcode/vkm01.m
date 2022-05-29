function [harmocs,residual]=vkm01(Data,freq,Fs)

    
%% 卡尔曼滤波器参数准备
r = 2000;     % 设置权重 2000
ford = 1;     % 设置滤波器阶数为1（也可以取2）
tol = 0.1;    % 误差精度
maxit = 500;  % 最大迭代次数
V_fp = kron(ones(size(Data,1),1),[freq freq*2 freq*3]); % 4 是叶片的个数, 5 是希望分离出的谐波个数
[nt,nch] = size(V_fp);


nt=size(Data,2);

parfor k=1:nt
    [xm,bwm,Tm,fl,rr,it,rv,xr] = vkm(Data(:,k),V_fp,Fs,r*ones(nch,1),ford,tol,maxit);
    residual(:,k) = Data(:,k) - sum(xr,2);
    harmocs(:,k)   = sum(xr,2);
end


end

filter