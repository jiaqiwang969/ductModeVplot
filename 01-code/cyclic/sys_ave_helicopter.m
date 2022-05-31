%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 级联循环维纳滤波器+卡尔曼滤波器
% 在噪声测量下对旋翼噪声进行提取，同时分离出谐波和宽带噪声
% 主要参考文献：
% (1) Liang Yu, Haijun Wu, Jerome Antoni, Weikang Jiang，
% Extraction and imaging of aerodynamically generated sound field of rotor 
% blades in the wind tunnel test, Mechanical Systems and Signal Processing, 
% Volume 116, 1 February 2019, Pages 1017-1028. 

% (2) Liang Yu, Jerome Antoni, Haijun Wu, Weikang Jiang, 
% Reconstruction of cyclostationary sound source based on a back-propagating
% cyclic wiener filter, Journal of Sound and Vibration,2019,442 :787-799.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
clear
close all

SampeFreq = 51200 ;  % 采样频率
% AnaChannelNums = 13;  % 分析第10通道13 24
SampleTime = 10; % 采样时间
Fs = SampeFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  循环维纳滤波器：对旋翼噪声进行提取（包含1阶和2阶统计信息）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 读取转速文件
% fid = fopen([DateFilePath,DateFileName,'-C-period.dat'], 'r');
% time_blade = fscanf(fid,'%f'); %传声器声压信号，每一圈的时间
% fclose(fid);
if 1
    [Datafilename,mypath1]=uigetfile('D:\A学习资料\程序与数据\绵阳数据3\Bo105\平飞\*.*','输入桨叶方位数据');
    S=strcat(mypath1,Datafilename);
    
    fid = fopen(S,'r','b');
    Nona=fread(fid,1,'single');
    Timenummic=fread(fid,1,'single');
    CLenthTime=fread(fid,1,'single');
    NLenthTime=fread(fid,1,'single');
    BladePositionSig1=fread(fid,CLenthTime*NLenthTime,'single');
    fclose(fid);
    
    
    Npiont=CLenthTime*NLenthTime/128;
%     Npiont=CLenthTime*NLenthTime
    for NNp=1:Npiont
        BladePositionSig(NNp)=sum(BladePositionSig1(((NNp-1)*128+1):NNp*128));
    end
    disp(['完成方位角数据读入!'])
    time_blade = BladePositionSig;
end 

% BO105 的信号读法
if 0
[Datafilename,mypath1]=uigetfile('D:\A学习资料\程序与数据\绵阳数据3\Bo105\平飞\*.*','输入桨叶方位数据'); 
S=strcat(mypath1,Datafilename);

fid = fopen(S,'r','b');
Nona=fread(fid,1,'single');
Timenummic=fread(fid,1,'single');
CLenthTime=fread(fid,1,'single');
NLenthTime=fread(fid,1,'single');
BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
fclose(fid);
disp(['完成方位角数据读入!'])
time_blade = BladePositionSig;
end


% 读取噪声数据
[Datafilename,mypath1]=uigetfile('D:\A学习资料\程序与数据\绵阳数据3\Bo105\平飞\*.*','输入试验数据'); 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  对原始数据简单时频分析
AnaChannelNums = 10;  % 分析第10通道13 24
y = RawData(1:SampleTime*SampeFreq,AnaChannelNums); % 从原始数据中读取分析数据
L = size(y,1);  % 分析数据的长度
Fs = SampeFreq;

% noise = y;
% signal = y = RawData(1:SampleTime*SampeFreq,1:10);
% SNR = 10*log10(var(signal)/var(noise))
% Lp = 20*log10(rms(signal)/(2*10^(-5)));

% 显示分析数据和基本的频谱分析(图1)
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
time_disp = 0.1;
frequency_disp = 6000;
% figure
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% set(findobj('type','axes'),'fontsize',28);
% xlabel('时间 (秒)'),ylabel('幅度'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('测量数据');

figure
subplot(211)
plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'k','LineWidth',0.5);
% xlabel('时间 (秒)'),ylabel('幅度'),
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);
% title('测量数据');
subplot(212)
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'k','LineWidth',1) % 只显示到6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
% title(['测量数据信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');


% 时频谱分析
Nw =  256;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
figure
spectrogram(y,Nw,Nv,Nfft, Fs,'yaxis')
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);

% 循环频谱相关性分析
alpha_max = 300; % 最大的循环频率 300
Nw = 256; % 窗口长度 256
opt.coh = 1;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);

figure,
plot(alpha(2:end),sum(abs(S(:,2:end))),'LineWidth',2),
% xlabel('循环频率 (Hz)','fontsize',12,'fontweight','b')
xlabel('Cyclic frequency (Hz)','fontsize',12,'fontweight','b')
ylabel('Sum of the spectral frequency (Hz)','fontsize',12,'fontweight','b');
% ylabel('谱频率求和 (Hz)','fontsize',12,'fontweight','b');
% title('循环频率显示');
set(findobj('type','axes'),'fontsize',12);
% xlim([0 alpha_max])

% % %  时域同步平均提取一阶循环平稳信号
% % % [~, syn_c_num] = size(time_blade);
% % Period_time = mean(time_blade);
% % Period_sample = round(Fs*Period_time*6);
% % Period_number = fix(L/Period_sample);
% % y_Period = y(1:Period_number*Period_sample);
% % y_Period_matrix = reshape(y_Period,Period_sample,Period_number);
% % My_Period = mean(y_Period_matrix,2);
% % ave_syn_y = repmat(My_Period,Period_number,1);
% % y_residual = y_Period - ave_syn_y;

% time_disp = 0.1;
% figure;
% subplot(311)
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('测量信号')
% subplot(312)
% plot([0:1/Fs:time_disp-1/Fs],ave_syn_y(1:time_disp*Fs,end),'LineWidth',0.5);
% % xlabel('时间 (秒)'),ylabel('幅度'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('提取的1阶循环平稳信号')
% subplot(313)
% plot([0:1/Fs:time_disp-1/Fs],y_residual(1:time_disp*Fs,end),'LineWidth',0.5);
% % xlabel('时间 (秒)'),ylabel('幅度'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('残差信号')
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');

time_disp = 0.5;
figure
plot([0:1/Fs:time_disp-1/Fs],ave_syn_y(1:time_disp*Fs,end),'LineWidth',0.5);

Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
time_disp = 0.1;
frequency_disp = 5000;
figure
subplot(311)
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % 只显示到6000Hz
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['测量信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])

subplot(312)
[S,f] = pwelch(ave_syn_y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % 只显示到6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['提取的1阶循环平稳信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');

subplot(313)
[S,f] = pwelch(y_residual,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % 只显示到6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['残差信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');

y = y_residual;
L = length(y);
%%  循环维纳滤波器
Na = 10; % 循环频率的阶数
orders = (1:Na/2); % 循环频率的个数
orders = [orders;-orders]; 
orders = orders(:); % 总的循环频率个数

alpha = 34.2/Fs; %34.2/Fs;%136.8/Fs/4; % 第一阶循环频率 136.8/Fs
phi = 2*pi*alpha*(0:L-1)'; % 相位的索引
fi = alpha*ones(L,1); % 频率的索引

% 估计循环互谱矩阵
Nw = 2^8;
[Syya,Syaya] = CyclicSpecMat(y,phi,fi,orders,Nw);
f = Syaya.f*Fs;

% 构建循环维纳滤波器
Ns = 1; % 循环平稳声源的个数
[G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,1e-6);

% 基于Gabor变换进行滤波
[xhat,cx] = CyclicSpatialFiltering(y,G(1,:,:),phi,fi,orders); %%xhat是输出的宽带信号
[L_xhat,~] = size(xhat);
nhat = y(1:L_xhat) - xhat;  %%nhat是剩余的噪声信号

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 前一秒的测量数据与提取的数据进行对比
time_disp = 0.1;% 0.05;
figure;
% subplot(211);
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',2);
% % xlabel('时间 [秒])'),ylabel('幅度[Pa]');
% ylabel('Magnitude [Pa]','fontsize',12,'fontweight','b');
% xlabel('Time [Second]','fontsize',12,'fontweight','b');
% title('残差信号');
subplot(212);
plot([0:1/Fs:time_disp-1/Fs],xhat(1:time_disp*Fs,end),'LineWidth',2);
% xlabel('时间 [秒]'),ylabel('幅度[Pa]');
ylabel('Magnitude [Pa]','fontsize',12,'fontweight','b');
xlabel('Time [Second]','fontsize',12,'fontweight','b');
title('进一步提取的二阶循环平稳信号');
% subplot(313);
% plot([0:1/Fs:time_disp-1/Fs],nhat(1:time_disp*Fs,end),'LineWidth',2);
% % xlabel('时间 (秒)'),ylabel('幅度');
% ylabel('Magnitude [Pa]','fontsize',12,'fontweight','b');
% xlabel('Time [Second]','fontsize',12,'fontweight','b');
% title('残差信号时间域表达');
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');

% % FFT psd
% frequency_disp = Fs/2;  %显示频率到**Hz
% f = Fs*(0:(L/2))/L;
% Y = fft(xhat); %FFT
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

% Nw =  Fs;   % window
% Nv =  ceil(3/4*Nw); % overlap
% Nfft = Nw;  % nfft
% frequency_disp = 6000;
% figure
% [S,f] = pwelch(residual,Nw,Nv,Nfft,Fs);
% plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'k','LineWidth',1.25) 
% ylabel('PSD [dB/Hz]','fontsize',12);
% xlabel('Frequency [Hz]','fontsize',12);
% title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
% xlim([0 frequency_disp])
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% set(gcf,'position',[400 250 500 250])

% figure
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% set(findobj('type','axes'),'fontsize',28);
% xlabel('时间 (秒)'),ylabel('幅度'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('测量数据');
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
time_disp = 0.1;
frequency_disp = 6000;
plot([0:1/Fs:time_disp-1/Fs],xhat(1:time_disp*Fs,end),'k','LineWidth',0.5);
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);
figure
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % 只显示到6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
% title(['测量数据信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

%% 数据分析
% 滤波后的循环频谱相关性分析
alpha_max = 300; % 最大的循环频率
Nw = 2^8; % 窗口长度
opt.coh = 1;

figure
[S,alpha,f,STFT,t,Nv] = Fast_SC(xhat(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end))),'k','LineWidth',2), hold on;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end))),'r--','LineWidth',1);
% [S,alpha,f,STFT,t,Nv] = Fast_SC(nhat(:,end),Nw,alpha_max,Fs,opt);
% plot(alpha(2:end),sum(abs(S(:,2:end))),'y','LineWidth',2);
% xlabel('循环频率 (Hz)')
xlabel('Cyclic frequency (Hz)','fontsize',12,'fontweight','b')
ylabel('Sum of the spectral fruequcy (Hz)','fontsize',12,'fontweight','b');
% legend('提取的信号','测量的信号','残差信号')
legend({'Extracted signal','Measured signal'},'FontSize',12)
% title('提取信号的循环频率显示');
set(findobj('type','axes'),'fontsize',12);

% 测量信号和提取信号的频谱分析
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
frequency_disp = 800;
figure
[S,f] = pwelch(xhat(1:2*Fs),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',2),hold on
[S,f] = pwelch(y(1:2*Fs),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'r','LineWidth',2),hold on
xlabel('频率 (Hz)'),ylabel('dB')
ylabel('PSD (dB/Hz)','fontsize',12,'fontweight','b');
xlabel('Time [Second]','fontsize',12,'fontweight','b');
legend('提取的数据','测量数据')
title('信号的频谱分析');
set(findobj('type','axes'),'fontsize',12);

% 时频谱分析
Nw =  256;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
figure
spectrogram(xhat,Nw,Nv,Nfft, Fs,'yaxis')
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  卡尔曼滤波器：分离出谐波和宽带噪声
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 读取转速文件
% fid = fopen([DateFilePath,DateFileName,'-C-period.dat'], 'r');
% time_blade = fscanf(fid,'%f'); %传声器声压信号，每一圈的时间
% fclose(fid);
if 1
    [Datafilename,mypath1]=uigetfile('E:\reseach\绵阳数据3\*.*','输入桨叶方位数据');
    S=strcat(mypath1,Datafilename);
    
    fid = fopen(S,'r','b');
    Nona=fread(fid,1,'single');
    Timenummic=fread(fid,1,'single');
    CLenthTime=fread(fid,1,'single');
    NLenthTime=fread(fid,1,'single');
    BladePositionSig1=fread(fid,CLenthTime*NLenthTime,'single');
    fclose(fid);
    
    
    Npiont=CLenthTime*NLenthTime/128;
%     Npiont=CLenthTime*NLenthTime
    for NNp=1:Npiont
        BladePositionSig(NNp)=sum(BladePositionSig1(((NNp-1)*128+1):NNp*128));
    end
    disp(['完成方位角数据读入!'])
    time_blade = BladePositionSig;
end 

% BO105 的信号读法
if 0
[Datafilename,mypath1]=uigetfile('D:\绵阳数据3\*.*','输入桨叶方位数据'); 
S=strcat(mypath1,Datafilename);

fid = fopen(S,'r','b');
Nona=fread(fid,1,'single');
Timenummic=fread(fid,1,'single');
CLenthTime=fread(fid,1,'single');
NLenthTime=fread(fid,1,'single');
BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
fclose(fid);
disp(['完成方位角数据读入!'])
time_blade = BladePositionSig;
end

% 根据转速信息同步一阶瞬时频率
f_rot=[]; % 瞬时频率
for i=1:length(time_blade)
f_i=1/time_blade(i); % 轴频，1/time_blade(i) 是叶频，4是桨叶的个数
f_rot=[f_rot f_i*ones(1,round(time_blade(i)*Fs))];
end
f_rot = f_rot.';

% 信号准备
ratio = 2;
y_xhat = xhat(1:ratio*Fs); 
% y_xhat = xhat(1:2*Fs);   %取部分信号送进卡尔曼滤波器进行滤波
fp = f_rot(1:size(y_xhat,1)); 
% t = (1:size(y_xhat,2))/Fs;                                 %离散时间信号

%% 卡尔曼滤波器参数准备
r = 2000;   % 设置权重 2000
ford = 1;    % 设置滤波器阶数为1（也可以取2）
tol = 0.1;    % 误差精度
maxit = 500;  % 最大迭代次数

V_fp = [4*fp 2*4*fp 3*4*fp 4*4*fp 5*4*fp 6*4*fp]; % 4 是叶片的个数, 5 是希望分离出的谐波个数
[nt,nch] = size(V_fp);

[xm,bwm,Tm,fl,rr,it,rv,xr] = vkm(y_xhat,V_fp,Fs,r*ones(nch,1),ford,tol,maxit);
figure,plot(0:maxit, rv,'LineWidth',2);
ylabel('residual error','fontsize',12,'fontweight','b');
xlabel('Iteration','fontsize',12,'fontweight','b');
%% 提取的谐波信号的表达
figure;
subplot(511)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,1),'LineWidth',2); 
legend({'叶频第1谐波'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(512)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,2),'LineWidth',2)
legend({'叶频第2谐波'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(513)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,3),'LineWidth',2)
legend({'叶频第3谐波'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(514)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,4),'LineWidth',2)
legend({'叶频第4谐波'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(515)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,5),'LineWidth',2)
legend({'叶频第5谐波'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');

%% 从旋翼噪声中恢复谐波信号（频谱）
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
frequency_disp = 1500;
figure
[S,f] = pwelch(y_xhat,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,2),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,3),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,4),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,6),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,5),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
legend({'宽频+谐波','叶频第1谐波','叶频第2谐波','叶频第3谐波','叶频第4谐波','叶频第5谐波','叶频第6谐波'},'FontSize',12)
%%

% 从旋翼噪声中恢复谐波信号（频谱）
% frequency_disp = 800;
% figure
% Spec_y = fft(y(1:ratio*Fs,1));
% f = linspace(0,Fs,length(Spec_y));
% plot(f,10*log10(abs(Spec_y/size(Spec_y,2))));
% xlim([0 frequency_disp])
% hold on
% Spec_xr1 = fft(xr(:,1));
% f = linspace(0,Fs,length(Spec_y));
% plot(f,10*log10(abs(Spec_xr1/size(Spec_xr1,2))));
% xlim([0 frequency_disp])
% Spec_xr2 = fft(xr(:,2));
% plot(f,10*log10(abs(Spec_xr2/size(Spec_xr2,2))));
% xlim([0 frequency_disp])
% Spec_xr3 = fft(xr(:,3));
% plot(f,10*log10(abs(Spec_xr3/size(Spec_xr3,2))));
% xlim([0 frequency_disp])
% Spec_xr4 = fft(xr(:,4));
% plot(f,10*log10(abs(Spec_xr4/size(Spec_xr4,2))));
% xlim([0 frequency_disp])
% Spec_xr5 = fft(xr(:,5));
% plot(f,10*log10(abs(Spec_xr5/size(Spec_xr5,2))));
% xlim([0 frequency_disp])


% 分离出来的谐波信号和宽带信号
residual = y_xhat - sum(xr,2);
harmocs  = sum(xr,2);
figure,
subplot(211)
plot(0:1/Fs:0.1-1/Fs,harmocs(1:0.1*Fs))
legend({'叶频谐波噪声'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(212)
plot(0:1/Fs:0.1-1/Fs,xhat(1:0.1*Fs),'LineWidth',0.5);
% legend({'宽频噪声'},'FontSize',12);
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);

plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% xlabel('时间 (秒)'),ylabel('幅度'),
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);
% % (3) 功率谱分析（Welch’s power spectral density estimate）
% Nw =  Fs/2;   % window
% Nv =  ceil(3/4*Nw); % overlap
% Nfft = Nw;  % nfft
% frequency_disp = Fs/4;
% figure
% [S,f] = pwelch(residual,Nw,Nv,Nfft,Fs);
% plot(f(1:frequency_disp/2),10*log10(S(1:frequency_disp2/)),'LineWidth',1) % 只显示到6000Hz
% % xlabel('Frequency (Hz)'),ylabel('dB'),
% ylabel('PSD [dB/Hz]','fontsize',12);
% xlabel('Frequency [Hz]','fontsize',12);
% % title(['提取的1阶循环平稳信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% (3) 功率谱分析（Welch’s power spectral density estimate）
Nw =  Fs/2;   % window
Nv =  ceil(3/4*Nw); % overlap
Nfft = Nw;  % nfft
frequency_disp = Fs/8;
figure
[S,f] = pwelch(xhat,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp/2),10*log10(S(1:frequency_disp/2)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
% title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
ylim([-80 0])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 200])

% 测量信号、谐波信号和宽带信号的频谱
frequency_disp = 2000;
figure
[S,f] = pwelch(y_xhat,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'b','LineWidth',2)
hold on
[S,f] = pwelch(harmocs,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'r--','LineWidth',2)
hold on
[S,f] = pwelch(residual,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'y--','LineWidth',2)
legend({'测量信号','叶频谐波噪声','宽频噪声'},'FontSize',12)
set(findobj('type','axes'),'fontsize',12);
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');

