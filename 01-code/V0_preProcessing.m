% Aim: use data to plot Vtu
% 2022-05-15 
% Ref: Behn M, Pardowitz B, Tapken U. Separation of tonal and broadband 
% noise components by cyclostationary analysis of the modal sound field in 
% a low-speed fan test rig[C]//International Conference of Fan Noise, 
% Aerodynamics, Applications and Systems. 2018: 18-20.
% wjq - 2022-05-17


clc;
clear;
close all;

%% add subfunction
addpath(genpath('.'));
chemin = '../database/01-rotateMic';

%% add Basic parameters

zH = 0.4;         % 测试距离
nk = 12;          % 传声器的数量
NumSM= 30;        % 测量的次数
a=0.185;          % 管道半径
S=pi*a^2;         % 管口面积
Fs = 102400 ;     % 采样频率
time=5;           % 采样时间

%% data processing
L_signal = Fs*time;             %信号长度
L_seg = round(L_signal/100);    %确定对信号处理划窗长度
Wind = hamming(L_seg);          %确定对数据进行汉宁窗处理
Noverlap = round(L_seg/2);      %确定信号划窗重叠率
Nfft = 2^(ceil(log2(L_seg))+1); %确定分析频率

rotor_speed=12000;               %轴转速信息
nk=12;

%% data processing
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
for i_file =Num_file
    eval(['load ''',chemin,'/','RotaryTest-12000-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
    Tdata{i_file}=Data(:,1:13);
    % Step01: 通过key signal将其分段,生成多个block，每个block 6 round, 历时 25 秒。
    % 在这里需要增加等角度采样的操作：
    [key_pulse,rotor_speed]=keyRotation(Data(:,14),Fs);
    cut_number(i_file)=floor(length(key_pulse)/nk)-1;
    data_resample_interval(i_file)=key_pulse(nk+1)-key_pulse((1));  
    for kb=1:cut_number(1)
        tmp=Tdata{i_file}(key_pulse((1+(kb-1)*nk)):key_pulse(1+(kb*nk)),:);
        data_block{kb,i_file}=resample(tmp,data_resample_interval(1),size(tmp,1));
    end
end
cut_number=cut_number(1);
    % Step02: ensember average 得到tonal noise, 历时 2.707163 秒。
      data_block_3d = reshape(cell2mat(data_block.'),data_resample_interval(1)*NumSM,13,cut_number);
      data_tonal_rms=mean(data_block_3d,3);
      data_tonal_rms2=mat2cell(data_tonal_rms,data_resample_interval(1)*ones(NumSM,1),[13]).'; % 形式与Tdata保持一致
    % Step03: r(t)=p(t)-s(t)
      data_tonal=kron(ones(cut_number,1),cell2mat(data_tonal_rms2));
      data_broadband=cell2mat(data_block)-data_tonal;



%% 作图

[Gx0,Gxx0,Fx0] = avgGxx('hann',50,'ACF',10,Fs,3200,Tdata{1, 1}(:,1)); %暂时fs手动微调
[Gx1,Gxx1,Fx1] = avgGxx('hann',50,'ACF',10,Fs,3200,data_tonal(:,1)); 
[Gx2,Gxx2,Fx2] = avgGxx('hann',50,'ACF',10,Fs,3200,data_broadband(:,1)); 
abs_q0=20*log10(abs(Gx0)/(2*10-5));
abs_q1=20*log10(abs(Gx1)/(2*10-5));
abs_q2=20*log10(abs(Gx2)/(2*10-5));
figure;plot(Fx0,abs_q0,'k','LineWidth',5);hold on;plot(Fx1,abs_q1,'b','LineWidth',2);   plot(Fx2,abs_q2,'r','LineWidth',2); 
xlim([0 50000])
grid on
grid minor
xlabel('Frequency/Hz')
ylabel('Sound pressure level/dB')




%% 分析数据
y = Tdata{1, 1}(:,1);%RawData(1:SampleTime*SampeFreq,AnaChannelNums); % 从原始数据中读取分析数据
L = size(y,1);  % 分析数据的长度
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
fq=rotor_speed/60*29*2.5;
plot(f(1:fq),10*log10(S(1:fq))) % 只显示到6000Hz
xlabel('频率 (赫兹)'),ylabel('dB'),
title(['信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])

% 循环频谱相关性分析
alpha_max = fq; % 最大的循环频率
Nw = 2^8; % 窗口长度
opt.coh = 1;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);
figure
plot(alpha(2:end),sum(abs(S(:,2:end)))),
xlabel('循环频率 (Hz)')
title('循环频率显示');

ave_syn_y = data_tonal(:,1);
y_residual = data_broadband(:,1);

%%  卡尔曼滤波器（去除一阶谐波信号）
% 功率谱分析（Welch’s power spectral density estimate）
Nw =Fs/20;   % window
Nv =  ceil(3/4*Nw); % overlap
Nfft = Nw;  % nfft
frequency_disp = Fs/2;
figure
[S,f] = pwelch(data_broadband(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp/20),20*log10(S(1:frequency_disp/20)/(2*10^-5)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])
hold on
[S,f] = pwelch(data_tonal(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp/20),20*log10(S(1:frequency_disp/20)/(2*10^-5)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])
legend('broadband','tonal');

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

% 滤波后的循环频谱相关性分析
alpha_max = fq; % 最大的循环频率
Nw = 2^8;       % 窗口长度
opt.coh = 1;



figure,
[S,alpha,f,STFT,t,Nv] = Fast_SC(xhat(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end)))), hold on;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y_residual(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end))),'r');
xlabel('循环频率 (Hz)');
legend('提取的信号','测量的信号');
title('提取信号的循环频率显示');

%% 测量信号和提取信号的频谱分析
Nw =   Fs;
Nv =   ceil(3/4*Nw);
Nfft = Nw;
figure
[S,f] = pwelch(xhat(:,end),Nw,Nv,Nfft,Fs);
plot(f(1:fq),20*log10(S(1:fq)/(2*10^-5)));hold on;
[S,f] = pwelch(y(:,end),Nw,Nv,Nfft,Fs);
plot(f(1:fq),20*log10(S(1:fq)/(2*10^-5)))
[S,f] = pwelch(data_broadband(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:fq),20*log10(S(1:fq)/(2*10^-5)))
[S,f] = pwelch(data_tonal(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:fq),20*log10(S(1:fq)/(2*10^-5)))
xlabel('频率 (Hz)'),ylabel('dB')
legend('提取的数据','测量数据','broadband','tonal')
title('信号的频谱分析');

