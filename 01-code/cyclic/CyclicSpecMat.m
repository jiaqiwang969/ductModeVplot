function [Syya,Syaya] = CyclicSpecMat(y,phi,alpha,orders,Window,Noverlap)
% [Syya,Syaya] = CyclicSpecMat(y,phi,alpha,orders,Window,Noverlap)
% 估计循环互功率谱矩阵Syya 和 Syaya
% 其中 y 是原始的测量信号，ya 是通过对循环频率alpha扩展的扩展信号

% 输入参数:
% y 是输入信号的矩阵，其中每一列包含一个通道的数据
% phi = 瞬时相位，单位是弧度
% alpha = 瞬时频率 (单位为Hz)，并通过除以采样频率归一化 = (d phi(t)/dt)/(2 pi*Fs)
% orders = 循环频率的阶数
% 注意: 如果这里输入的 phi 和 alpha 是矩阵, 说明每一列对应不同的循环相位和循环频率;
% 否则, alpha 的阶数直接从输入向量 'orders'中选取
% Window = Welch 方法中使用的窗口长度
% 如果窗口长度的输入是标量, 则表明为汉宁窗的长度
% Noverlap (optional) = 按照理论一般选取为 3/4*Nw, 即3/4的窗口长度

% 输出参数:
% Syya = 循环互功率谱矩阵, (y 和 ya 之间的互谱矩阵)
% Syaya = 循环互功率谱矩阵, (ya 和 ya 之间的互谱矩阵)
% 其中 Syya 和 Syaya 是结构类型的数据:
%                   Syya.S = 循环互功率谱矩阵
%                   Syya.f = 频率向量
%                   Syya.K = blocks块的个数
%                   Syya.Var_Reduc = 方差减少因子
% 
% 该函数调用 'CSM_W' 子函数.

[L,M] = size(y);

if length(Window) == 1
    Nw = Window;
    Window = hanning(Nw);
else
    Nw = length(Window);
end
nfft = Nw;
if nargin < 6
    Noverlap = fix(3/4*Nw);
end


if min(size(phi)) > 1
    Harmonics = 0;
    Na = size(phi,2);
else
    Harmonics = 1;
    Na = length(orders);
end

% ya = [];
% for ka = 1:Na
%     ya = [ya y.*repmat(exp(orders(ka)*1i*phi(:)).*alpha(:),1,M)];
% end
ya = zeros(L,Na*M);
for ka = 1:Na
    if Harmonics == 1
         temp = exp(orders(ka)*1i*phi(:)).*alpha(:);
        % temp = exp(orders(ka)*1i*phi(:));
    else
        temp = exp(1i*phi(:,ka)).*alpha(:,ka);
    end
    for m = 1:M
        ya(:,(ka-1)*M+m) = y(:,m).*temp;
    end
end

Syya = CSM_W(y,ya,nfft,Noverlap,Window);
Syaya = CSM_W(ya,ya,nfft,Noverlap,Window);