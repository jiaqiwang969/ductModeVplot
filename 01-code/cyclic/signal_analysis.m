%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%信号分析%%%%%%%%%%%%%%%
% 主要分析方法：
% 1.基础分析方法
% (1) 信号时域波形
% (2) FFT 频谱分析
% (3) PSD 功率谱分析
% (4) 包络谱分析
% (5) STFT 时频分析
% (6) Spectral Kurtosis/Kurtogram 谱峭度
% 2.循环平稳分析
% (7) SC 循环谱相关/循环谱相干分析(Fast-SC算法)
% 3.模态分析
% (8) EMD 经验模态分解
%%%%%%%%%%%%欢迎继续补充%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear,close all;
SampeFreq =51200 ;  % 采样频率(Hz)
% SampleTime = 30; % 采样时间(s)
Fs = SampeFreq;
%% 读取分析数据
[DataFileName,DataFilePath]=uigetfile('','选择原始数据'); % 打开读取数据窗口
load(strcat(DataFilePath,DataFileName));
disp('完成原始数据读入!')
RawData =xhat;    % x为导入数据中的变量名，将其放入RawData中
[m,n] = size(RawData);
if m<n
    RawData = RawData'; % 将RawData统一改写为 (采样点数*通道数) 的形式
end
AnaChannelNums = 30;  % 分析的通道数
AnaTime = 1;   % 分析的时间长度(s)
y = RawData(1:AnaTime*SampeFreq,AnaChannelNums); % 从原始数据中读取分析数据
L = size(y,1);  % 分析数据的长度
%% 基础分析方法
% (1) 分析数据的时域波形
% figure
% t = 0:1/Fs:AnaTime-1/Fs;    % 定义时间变量
% plot(t,y,'k','LineWidth',1.5);
% time_disp = 10;    % 显示时间长度
% xlim([0 time_disp])
% ylabel('Amplitude');    % 注意：需加上单位[m/s^2]或[Pa]
% xlabel('Time [Second]');
% title('Time Domain Waveform','fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% set(gcf,'position',[400 250 500 250])
% % 
% % (2) FFT频谱分析
frequency_disp = Fs/2;  %显示频率到**Hz
f = Fs*(0:(L/2))/L;
Y = fft(y); %FFT
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure
plot(f,P1,'k','LineWidth',1.5) 
title('FFT Amplitude Spectrum','fontweight','b')
xlabel('f [Hz]')
ylabel('Amplitude')    % 注意：需加上单位[m/s^2]或[Pa]
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])


% (3) 功率谱分析（Welch’s power spectral density estimate）
Nw =  Fs;   % window
Nv =  ceil(3/4*Nw); % overlap
Nfft = Nw;  % nfft
frequency_disp = Fs/2;
figure
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

% (4) 包络谱分析
frequency_disp = 2e3;  %显示频率到2000Hz
f = Fs*(0:(L/2))/L;
y1=y-mean(y);   %去除直流分量
Yh=abs(hilbert(y1));    %希尔伯特变换
lnY=length(Yh);
Y1=fft(Yh-mean(Yh))/lnY;    %包络的FFT
Y1=abs(Y1(1:L/2+1));
figure
plot(f,Y1,'k','LineWidth',1.5);
title('Envelop Spectrum','fontweight','b')
xlim([0 frequency_disp]) %频率的显示范围
xlabel('f [Hz]')
ylabel('Amplitude')    % 注意：需加上单位[m/s^2]或[Pa]
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

% 
% % (5)STFT时频谱分析（Spectrogram using short-time Fourier transform）
frequency_disp = Fs/2;    % 显示频率长度
Nw =  2^11;  % 窗长
window = hanning(Nw);      % hanning窗
Nv =  ceil(3/4*Nw); % overlap:重叠样本数:3/4*Nw
Nfft = 4*Nw;  % nfft:DFT点数
[S, F, T] = spectrogram(y, window, Nv, Nfft, Fs);
figure
imagesc(T, F, log10(abs(S)))
set(gca, 'YDir', 'normal')
xlabel('Time [Seconds]')
ylabel('Frequency [Hz]')
title('STFT Spectrum','fontweight','b')
h=colorbar;
set(get(h,'Title'),'string','dB/Hz');
ylim([0 frequency_disp]) %频率的显示范围
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);
set(gcf,'position',[400 250 500 350])
% 
% % (6) Spectral Kurtosis/Kurtogram 谱峭度
% figure
% pkurtosis(y,Fs) % Spectral Kurtosis
% title('Spectral Kurtosis','fontweight','b')
% ylabel('Spectral Kurtosis'),xlabel('f [kHz]'),
% xlim([0 Fs/2000])
% set(findobj('type','axes'),'fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(gcf,'position',[400 250 500 250])
% 
% figure
% kurtogram(y,Fs) % Kurtogram
% xlim([0 Fs/2000])
% set(findobj('type','axes'),'fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(gcf,'position',[400 250 500 400])
%% 循环平稳分析
% (7) SC 循环谱相关/循环谱相干分析  (Fast-SC算法)
% 运行需调用子函数 Fast_SC.m
alpha_max = 800; % 最大的循环频率
Nw = 2^8; % 窗口长度
opt.coh = 1;    %为0时计算Spectral Correlation，为1时计算Spectral Coherence
[S,alpha,f,STFT,~,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);  %利用Fast-SC算法计算SC

% 绘制谱相关/谱相干
figure
surf(alpha(2:end),f,abs(S(:,2:end))),axis xy,colorbar()
shading flat
if opt.coh == 0,title('Spectral Correlation','fontweight','b'),else,title('Spectral Coherence','fontweight','b'),end
xlabel('cyclic frequency \alpha [Hz]'),ylabel('f [Hz]'),
xlim([0 alpha_max]),ylim([0 Fs/2])
view(30,30)
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);
set(gcf,'position',[400 250 500 350])
% 绘制平均谱相关/平均谱相干
figure
plot(alpha(2:end),mean(abs(S(:,2:end))),'k'),xlim([0 alpha_max])
if opt.coh == 0,title('Mean Spectral Correlation','fontweight','b'),else,title('Mean Spectral Coherence','fontweight','b'),end
xlabel('cyclic frequency \alpha (Hz)')
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);
set(gcf,'position',[400 250 500 250])
% %% 模态分析
% % (8) EMD 经验模态分解
% [imf_emd,residual]=emd(y,'MaxNumIMF',5);
% imf_emd = imf_emd';
% figure;
% frequency_disp = Fs/2;  %显示频率到**Hz
% subplot(size(imf_emd,1)+1,2,1);
% t = 0:1/Fs:AnaTime-1/Fs;    % 定义时间变量
% plot(t,y,'k','LineWidth',1.5);
% title('EMD');
% ylabel('Signal')
% subplot(size(imf_emd,1)+2,2,2);
% f = Fs*(0:(L/2))/L;
% Y = fft(y); %FFT
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% plot(f,P1,'k','LineWidth',1.5);
% xlim([0 frequency_disp]) %频率的显示范围
% title('Amplitude Spectrum');
% ylabel('Signal')
% for i = 2:size(imf_emd,1)+1
%     subplot(size(imf_emd,1)+2,2,i*2-1);
%     plot(t,imf_emd(i-1,:),'k','LineWidth',1.5);
%     ylabel(['IMF',num2str(i-1)])
%     subplot(size(imf_emd,1)+2,2,i*2);
%     Y = fft(imf_emd(i-1,:)); %FFT
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     plot(f,P1,'k','LineWidth',1.5);
%     ylabel(['IMF',num2str(i-1)])
%     xlim([0 frequency_disp]) %频率的显示范围
% end
% subplot(size(imf_emd,1)+2,2,2*size(imf_emd,1)+3);
% plot(t,residual,'k','LineWidth',1.5);
% xlabel('Time [s]')
% ylabel('Residual')
% subplot(size(imf_emd,1)+2,2,2*size(imf_emd,1)+4);
% Y = fft(residual); %FFT
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% plot(f,P1,'k','LineWidth',1.5);
% xlabel('Frequency [Hz]')
% ylabel('Residual')
% xlim([0 frequency_disp]) %频率的显示范围
% set(findobj('type','axes'),'fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(gcf,'position',[300 -400 800 850])