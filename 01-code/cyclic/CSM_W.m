function Mat = CSM_W(y,x,nfft,Noverlap,Window)     
% Mat = CSM_W(y,x,nfft,Noverlap,Window)
% 信号 x 和 信号 y 之间的互谱矩阵估计(矩阵)
%
% 输入:
% y, x = 信号矩阵, 其中每一列是一个通道的信号
% Window = 窗口长度
% nfft = 每个窗口做FFT点的个数 推荐：nfft = 2*NWind, 即窗口长度NWind的两倍
% Noverlap = 窗口重叠 推荐：Noverlap = 1/2*Nwind (半正弦)， Noverlap = 2/3*Nwind (汉宁窗)
% CSM_W 调用 'CPS_W' 函数.
%
% 输出: 
% S 是一个结构类型的数据:
%                   Mat.S = 互谱矩阵, 其中维度为 (x列)x(y列)x(Nwind/2+1)
%                   Mat.f  = 频率向量
%                   Mat.K = 快拍block的个数
%                   Mat.Var_Reduc = 方差减小因子


Kx = size(x,2);
Ky = size(y,2);

Nw = length(Window);
if Nw == 1
    Nw = Window;
end
Mat.S = zeros(Ky,Kx,Nw/2+1);
count = 1;

same = 0;
if Kx == Ky
    if all(all(y==x))
        same = 1;
    end
end
if same == 1
    for i = 1:Kx
        Sij = CPS_W(y(:,i),x(:,i),nfft,Noverlap,Window);
        Mat.S(i,i,:) = real(Sij.S(1:Nw/2+1));
        for j = 1:i-1
            waitbar(count/Kx/(Kx-1)*2),count = count + 1;
            Sij = CPS_W(y(:,i),x(:,j),nfft,Noverlap,Window);
            Mat.S(i,j,:) = Sij.S(1:Nw/2+1);
            Mat.S(j,i,:) = conj(Sij.S(1:Nw/2+1));
        end
    end
else
    for i = 1:Ky
        for j = 1:Kx
            waitbar(count/Kx/Ky),count = count + 1;
            Sij = CPS_W(y(:,i),x(:,j),nfft,Noverlap,Window);
            Mat.S(i,j,:) = Sij.S(1:Nw/2+1);
        end
    end
end

Mat.Var_Reduc = Sij.Var_Reduc; 
Mat.K = Sij.K;
Mat.f = Sij.f; 


function Spec = CPS_W(y,x,nfft,Noverlap,Window)
% Spec = CPS_W(y,x,Noverlap,Window)
% 信号 x 和 y 之间的互谱估计 (向量)
% 输入:
% y, x = 输入信号
% Window = 窗口长度
% nfft = 每个窗口做FFT点的个数 推荐：nfft = 2*NWind, 即窗口长度NWind的两倍
% Noverlap = 窗口重叠 推荐：Noverlap = 1/2*Nwind (半正弦)， Noverlap = 2/3*Nwind (汉宁窗)
%
% 输出:
% Spec 是结构类型的数据输出:
%                   Spec.S = 互谱向量
%                   Spec.f =  频率(向量)
%                   Spec.K = 快拍block的个数
%                   Spec.Var_Reduc = 方差减小因子

if length(Window) == 1
    Window = hanning(Window);
end
Window = Window(:);
n = length(x);          % 数据点个数
nwind = length(Window); % 窗口长度

% 输入检查
if nwind <= Noverlap,error('Window length must be > Noverlap');end
if nfft < nwind,error('Window length must be <= nfft');end

y = y(:);
x = x(:);
K = fix((n-Noverlap)/(nwind-Noverlap));	% 窗口个数

% 计算 CPS 谱
index = 1:nwind;
f = (0:nfft-1)/nfft;
CPS = 0;

for i=1:K
    Yw = fft(Window.*y(index),nfft);		% Yw(f)
    Xw = fft(Window.*x(index),nfft);		% Xw(f-a)
    CPS = Yw.*conj(Xw) + CPS;
    index = index + (nwind - Noverlap);
end

% 归一化
KMU = K*norm(Window)^2;	% 归一化因子, 保证渐近无偏
CPS = CPS/KMU;

% 计算方差减小因子
Window = Window(:)/norm(Window);
Delta = nwind - Noverlap;
R2w = xcorr(Window);
k = nwind+Delta:Delta:min(2*nwind-1,nwind+Delta*(K-1));
if length(k) >1
    Var_Reduc = R2w(nwind)^2/K + 2/K*(1-(1:length(k))/K)*(R2w(k).^2);
else
    Var_Reduc = R2w(nwind)^2/K;
end

% 输出参数
Spec.S = CPS;
Spec.f = f;
Spec.K = K;
Spec.Var_Reduc = Var_Reduc;





