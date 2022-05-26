function [c,E] = rgt(x,g,R)
% c = rgt(x,g,R)
% 对信号x进行Gabor 变换
% 输入： 
% x 是输入的信号；
% g为分析窗口长度；
% R为窗口函数每次移动的长度；
%输出： 
% c 为Gabor 变换后的系数矩阵，其中列的序号对应频率，行的序号对应时间
% E 为Gabor 变换后的系数矩阵与短时傅里叶变换后的系数矩阵之间的相位关系：c_Gabor = c_STFT.*E


x = x(:);
g = g(:);
L = length(x);
Nw = length(g);
if rem(Nw,2 > 0)
    error('the length of window g must be even !');
end
K = floor((L-Nw)/R+1);

c = zeros(Nw/2+1,K);
ind = 1;
wi = 2i*pi*(0:Nw/2)'/Nw;
for k = 1:K
    temp = fft(x(ind:ind+Nw-1).*conj(g));
    c(:,k) = temp(1:Nw/2+1).*exp(-wi*(k-1)*R);
    ind = ind + R;
end

if nargout > 1
    E = ones(Nw/2+1,K);
    for k = 2:K
        E(:,k) = exp(-wi*(k-1)*R);
    end
end