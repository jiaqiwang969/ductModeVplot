%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  循环维纳滤波器
%               在噪声测量下对（二阶循环平稳）旋翼噪声进行提取
% 主要参考文献：
% (1)Liang Yu, etc, Extraction and imaging of aerodynamically generated sound field of rotor blades 
% in the wind tunnel test, MSSP2018
% (2)Liang Yu, etc, Reconstruction of cyclostationary sound source based on a back-propagating
% cyclic wiener filter, JSV2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
clear
close all
% 文件路径
SampeFreq = 51200 ;  % 采样频率
AnaChannelNums = 1:135;  % 分析第10通道
SampleTime = 1; % 采样时间
% 读入原始测量数据
[Datafilename,mypath1]=uigetfile('D:\A学习资料\程序与数据\绵阳数据3\Bo105\悬停\*.*','输入桨叶麦克风测量数据'); 
S=strcat(mypath1,Datafilename);
fid = fopen(S,'r','b');
Binaryfs=fread(fid,1,'single');
Binarynummic=fread(fid,1,'single');
CLenth=fread(fid,1,'single');
NLenth=fread(fid,1,'single');
RawData=zeros(CLenth*NLenth,Binarynummic);
tic
for NNi=1:NLenth
    RawData(((NNi-1)*CLenth+1):(NNi*CLenth),:)=fread(fid,[CLenth Binarynummic ],'single');
end
toc
fclose(fid);
disp(['完成试验数据读入!'])
%% 分析数据
y = RawData(1:SampleTime*SampeFreq,AnaChannelNums); % 从原始数据中读取分析数据
L = size(y,1);  % 分析数据的长度
Fs = SampeFreq;
% 显示分析数据和基本的频谱分析
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
[S,f] = pwelch(y(:,end),Nw,Nv,Nfft,Fs);
figure
subplot(211)
plot([0:1/Fs:1-1/Fs],y(1:Fs,1));
xlabel('时间 (秒)'),ylabel('幅度'),
title('显示第一秒内的测量数据');
subplot(212)
plot(f(1:6000),10*log10(S(1:6000))) % 只显示到6000Hz
xlabel('频率 (赫兹)'),ylabel('dB'),
title(['信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])

% 循环频谱相关性分析
alpha_max = 300; % 最大的循环频率
Nw = 2^8; % 窗口长度
opt.coh = 1;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);
figure
plot(alpha(2:end),sum(abs(S(:,2:end)))),
xlabel('循环频率 (Hz)')
title('循环频率显示');

%%  时域同步平均（去除一阶谐波）
% [Datafilename,mypath1]=uigetfile('D:\A学习资料\程序与数据\绵阳数据3\Bo105\平飞\*.*','输入桨叶方位数据'); 
% S=strcat(mypath1,Datafilename);
% fid = fopen(S,'r','b');
% Nona=fread(fid,1,'single');
% Timenummic=fread(fid,1,'single');
% CLenthTime=fread(fid,1,'single');
% NLenthTime=fread(fid,1,'single');
% BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
% fclose(fid);
% disp(['完成方位角数据读入!'])
% time_blade = BladePositionSig;
% 
%  同步平均提取一阶循环平稳信号
% [~, syn_c_num] = size(time_blade);
% ave_syn_y = zeros(size(y));
% y_residual = zeros(size(y));
ave_syn_y = zeros(49698,135);
y_residual = zeros(49698,135);

for i = 1:135
    time_blade = 136;
    Period_time = 1/136;
    Period_sample = round(Fs*Period_time*6);
    Period_number = fix(L/Period_sample);
    y_Period = y(1:Period_number*Period_sample,i);
    y_Period_matrix = reshape(y_Period,Period_sample,Period_number);
    My_Period = mean(y_Period_matrix,2);
    ave_syn_y(:,i) = repmat(My_Period,Period_number,1);
    y_residual(:,i) = y_Period - ave_syn_y(:,i);
end
%%  卡尔曼滤波器（去除一阶谐波信号）
% % 读取转速文件
% [Datafilename,mypath1]=uigetfile('D:\A学习资料\程序与数据\绵阳数据3\Bo105\平飞\*.*','输入桨叶方位数据'); 
% S=strcat(mypath1,Datafilename);
% fid = fopen(S,'r','b');
% Nona=fread(fid,1,'single');
% Timenummic=fread(fid,1,'single');
% CLenthTime=fread(fid,1,'single');
% NLenthTime=fread(fid,1,'single');
% BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
% fclose(fid);
% disp(['完成方位角数据读入!'])
% time_blade = BladePositionSig;
%     
% % 根据转速信息同步一阶瞬时频率
% f_rot=[]; % 瞬时频率
% for i=1:length(time_blade)
%     f_i=1/time_blade(i); % 轴频，1/time_blade(i) 是叶频，4是桨叶的个数
%     f_rot=[f_rot f_i*ones(1,round(time_blade(i)*Fs))];
% end
% f_rot = f_rot.';

% % 信号准备
% f_rot =34;
% fp = f_rot*ones(size(y,1),1); 
% % t = (1:size(y_xhat,2))/Fs; %离散时间信号
% 
% % 卡尔曼滤波器参数准备
% ave_syn_y = zeros(size(y));
% y_residual = zeros(size(y));
% for i=1:135
%     r = 2000;   % 设置权重 2000
%     ford = 1;    % 设置滤波器阶数为1（也可以取2）
%     tol = 0.1;    % 误差精度
%     maxit = 500;  % 最大迭代次数
%     V_fp = [4*fp 2*4*fp 3*4*fp 4*4*fp 5*4*fp 6*4*fp 7*4*fp 8*4*fp 9*4*fp 10*4*fp 11*4*fp 12*4*fp 13*4*fp 14*4*fp 15*4*fp 16*4*fp 17*4*fp 18*4*fp 19*4*fp 20*4*fp 21*4*fp 22*4*fp 23*4*fp 24*4*fp 25*4*fp 27*4*fp 28*4*fp 29*4*fp 30*4*fp 31*4*fp 32*4*fp 33*4*fp 34*4*fp 35*4*fp 37*4*fp 38*4*fp 39*4*fp 40*4*fp  0*fp]; % 4 是叶片的个数, 5 是希望分离出的谐波个数
%     [nt,nch] = size(V_fp);
%     [xm,bwm,Tm,fl,rr,it,rv,xr] = vkm(y(:,i)',V_fp,Fs,r*ones(nch,1),ford,tol,maxit);  
%     ave_syn_y(:,i)=sum(xr(:,1:39),2);
%     y_residual(:,i) = y(:,i)- ave_syn_y(:,i);    
% end

%% 信号分析 
% % FFT频谱分析
% frequency_disp = Fs/2;  %显示频率到**Hz
% f = Fs*(0:(L/2))/L;
% Y = fft(y_residual(:,1)); %FFT
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% figure
% plot(f,P1,'k','LineWidth',1.5) 
% title('FFT Amplitude Spectrum','fontweight','b')
% xlabel('f [Hz]')
% ylabel('Amplitude')    % 注意：需加上单位[m/s^2]或[Pa]
% xlim([0 frequency_disp])
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% set(gcf,'position',[400 250 500 250])

% 功率谱分析（Welch’s power spectral density estimate）
Nw =Fs/20;   % window
Nv =  ceil(3/4*Nw); % overlap
Nfft = Nw;  % nfft
frequency_disp = Fs/2;
figure
[S,f] = pwelch(y_residual(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp/20),10*log10(S(1:frequency_disp/20)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

%%  循环维纳滤波（提取二阶循环平稳信号）
L1 = size(y_residual,1);
Na = 10; % 循环频率的阶数
orders = (1:Na/2); % 循环频率的个数
orders = [orders;-orders]; 
orders = orders(:); % 总的循环频率个数
alpha = 34.2/Fs; % 第一阶循环频率 /136.8
phi = 2*pi*alpha*(0:L1-1)'; % 相位的索引
fi = alpha*ones(L1,1); % 频率的索引

% 估计循环互谱矩阵
Nw = 2^8;
[Syya,Syaya] = CyclicSpecMat(y_residual(:,:),phi,fi,orders,Nw);
f = Syaya.f*Fs;

% 构建循环维纳滤波器
Ns = 1; % 循环平稳声源的个数
[G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,1e-6);

% 基于Gabor变换进行滤波
[xhat,cx] = CyclicSpatialFiltering(y_residual(:,:),G(:,:,:),phi,fi,orders);
% [L_xhat,~] = size(xhat);
% nhat = y_residual(1:L_xhat) - xhat;
% y_residual=y_residual(1:51136,:);
% nhat = y_residual - xhat;
%% 绘制图像
% 前一秒的测量数据与提取的数据进行对比
figure;
subplot(311);
plot([0:1/Fs:0.1-1/Fs],y(1:0.1*Fs,end));
xlabel('时间 (秒)'),ylabel('幅度'),
title('测量数据');
subplot(312);
plot([0:1/Fs:0.1-1/Fs],y_residual(1:0.1*Fs,end));
xlabel('时间 (秒)'),ylabel('幅度'),
title('去除谐波剩余信号');
subplot(313);
plot([0:1/Fs:0.1-1/Fs],xhat(1:0.1*Fs,end));
xlabel('时间 (秒)'),ylabel('幅度'),
title('提取的循环平稳信号');
% subplot(313);
% plot([0:1/Fs:0.1-1/Fs],nhat(1:0.1*Fs,end));
% xlabel('时间 (秒)'),ylabel('幅度'),
% title('残差信号时间域表达');

% 滤波后的循环频谱相关性分析
alpha_max = 300; % 最大的循环频率
Nw = 2^8; % 窗口长度
opt.coh = 1;

figure,
[S,alpha,f,STFT,t,Nv] = Fast_SC(xhat(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end)))), hold on;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y_residual(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end))),'r');
% [S,alpha,f,STFT,t,Nv] = Fast_SC(nhat(:,end),Nw,alpha_max,Fs,opt);
% plot(alpha(2:end),sum(abs(S(:,2:end))),'y');
xlabel('循环频率 (Hz)')
legend('提取的信号','测量的信号')
% legend('提取的信号','测量的信号','残差信号')
title('提取信号的循环频率显示');

% % 测量信号和提取信号的频谱分析
% Nw =  Fs;
% Nv =  ceil(3/4*Nw);
% Nfft = Nw;
% figure
% [S,f] = pwelch(xhat(:,end),Nw,Nv,Nfft,Fs);
% plot(f(1:5000),10*log10(S(1:5000))),hold on
% [S,f] = pwelch(y(:,end),Nw,Nv,Nfft,Fs);
% plot(f(1:5000),10*log10(S(1:5000)))
% % [S,f] = pwelch(nhat(:,1),Nw,Nv,Nfft,Fs);hold on
% % plot(f(1:5000),10*log10(S(1:5000)))
% xlabel('频率 (Hz)'),ylabel('dB')
% legend('提取的数据','测量数据')
% title('信号的频谱分析');

% figure
% [S,f] = pwelch(nhat(:,1),Nw,Nv,Nfft,Fs);hold on
% plot(f(1:5000),10*log10(S(1:5000)))
% title('信号的频谱分析');