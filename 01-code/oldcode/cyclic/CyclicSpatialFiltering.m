function [x,cx] = CyclicSpatialFiltering(y,G,phi,alpha,orders,Window,R)
% [x,cx] = CyclicSpatialFiltering(y,G,phi,alpha,orders,Window,R)
% 基于短时傅里叶变换的方式对循环平稳信号进行提取
% 
% 输入:
% y =  输入信号, 其中每一列为一个通道信号;
% G = 分离滤波器
% phi = 瞬时相位，单位是弧度
% alpha = 瞬时频率 (单位为Hz)，并通过除以采样频率归一化 = (d phi(t)/dt)/(2 pi*Fs)
% orders = 循环频率的阶数
% 注意: 如果这里输入的 phi 和 alpha 是矩阵, 说明每一列对应不同的循环相位和循环频率;
% 否则, alpha 的阶数直接从输入向量 'orders'中选取
% Window (Gabor 变换参数) 
% = Welch 方法中使用的窗口长度, 用来计算Gabor变换的系数,缺省为Hanning(Nw)
% R (可选,Gabor 变换参数) = 窗口移动长度, default = Nw/4, 对应重复的窗口为3/4*Nw
% 该函数调用 'CSM_W', 'rgt'  'irgt' 子函数.
%
% 输出:
% --------
% x = 滤波后的信号, 即提取的循环平稳成分
% cx = 循环平稳成分Gabor的系数

[L,M] = size(y);
N = size(G,1);

if nargin < 6
    Nw = size(G,3);
    Nw = 2*(Nw-1);
    R = fix(Nw/4);
    Window = hanning(Nw);
else
    Nw = length(Window);
end

if min(size(phi)) > 1
    Harmonics = 0;
    Na = size(phi,2);
else
    Harmonics = 1;
    Na = length(orders);
end

K = floor((L-Nw)/R+1);
cx = zeros(N,Nw/2+1,K);
for ka = 1:Na
    if Harmonics == 1
        temp = exp(orders(ka)*1i*phi(:)).*alpha(:);
    else
        temp = exp(1i*phi(:,ka)).*alpha(:,ka);
    end
    for m = 1:M
        cya = rgt(y(:,m).*temp,Window,R);
        for n = 1:N
            cx(n,:,:) = squeeze(cx(n,:,:)) + repmat(squeeze(G(n,(ka-1)*M+m,:)),1,K).*cya;
        end
    end
end

L = (K-1)*R + Nw;
x = zeros(L,N);
for n = 1:N
    x(:,n) = irgt(squeeze(cx(n,:,:)),Window,R);
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cya = zeros(Na*M,Nw/2+1,K);
for ka = 1:Na
    temp = exp(orders(ka)*1i*phi(:)).*alpha(:);
    for m = 1:M
        cya((ka-1)*M+m,:,:) = rgt(y(:,m).*temp,Window,R);
    end
end

cx = zeros(M,Nw/2+1,K);
for k = 1:Nw/2+1
    cx(:,k,:) = squeeze(G(:,:,k))*squeeze(cya(:,k,:));
end

L = (K-1)*R + Nw;
x = zeros(L,M);
for m = 1:M
    x(:,m) = irgt(squeeze(cx(m,:,:)),Window,R);
end


