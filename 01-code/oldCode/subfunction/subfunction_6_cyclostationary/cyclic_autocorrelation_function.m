function R = cyclic_autocorrelation_function(x,Nfft,max_tau)
%输入：x：一维时间序列
%      Nfft：用于快速傅立叶变换的数据点数，等于循环频率的数据点数
%      max_tau：最大时延点数
%输出：时间序列x的循环自相关函数矩阵R
%      R维数: [2* max_tau+1, Nfft]

[rows,cols] = size(x);
if rows > cols
  x = x';
end
x = x - mean(x);

n = floor((length(x)-2*max_tau-1)/Nfft);  % 确定求平均的次数
r = zeros(2*max_tau+1,Nfft);           % 时变自相关函数
temp = zeros(Nfft,n);
t = (1:Nfft*n); 
for k = -max_tau:max_tau   
   temp(:) = conj(x(t+max_tau)).*x(t+k+max_tau);  %顺序按列填满数据
   if n==1
        r(k+1+max_tau,:) = temp';
    else
        r(k+1+max_tau,:) = mean(temp');   %对每列数据求平均
   end 
end

R = zeros(2*max_tau+1,Nfft);
for k = -max_tau:max_tau
R(k+1+max_tau,:) = exp(-j*pi*((0:Nfft-1)/Nfft)*k).*fft(r(k+1+max_tau,:))/Nfft;  
%对称形式的循环自相关函数
end

R = fftshift(R,2);  %将循环频率以0为中心点对称放置
