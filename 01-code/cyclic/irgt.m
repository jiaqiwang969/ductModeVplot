function [x,gg] = irSTFT(c,g,R,I)
% [x,gg] = irSTFT(c,g,R,I)
% Gabor反变换
% 输入参数：
% c 是 Gabor变换产生的系数矩阵
% g 是 时域的窗口函数
% R 是窗口函数移动的距离（对应于正变换）
% 如果输入参数大于3，信号x只重建对应第 I 个频率点的信号
% 输出参数：gg 是信号边界处截断的数目

[Nw,K] = size(c);
Nw = 2*(Nw-1);
g = g(:);
g = g/sum(g(1:R:Nw).^2);
L = (K-1)*R + Nw;

if nargin < 4
    c = ifft([c;conj(c(Nw/2:-1:2,:))]);
else
    if isempty(I)
        c = ifft([c;conj(c(Nw/2:-1:2,:))]);
    else
        C = zeros(Nw,K);
        n = [0:Nw-1]';
        for k = I
            C = C + exp(2i*pi*(k-1)*n/Nw)*c(k,:)/Nw;
        end
        c = C;
        clear C
    end
end

x = zeros(L,1);
ind = 1:Nw;
for k = 1:K
    x(ind) = x(ind) + g.*c(mod(ind-1,Nw)+1,k);
    ind = ind + R;
end
x = real(x);

if nargout > 1
    gg = zeros(L,1);
    ind = 1:Nw;
    for k = 1:K
        gg(ind) = gg(ind) + g.^2;
        ind = ind + R;
    end
end
