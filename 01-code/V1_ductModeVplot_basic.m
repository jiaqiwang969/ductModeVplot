% 8000-rmp——vplot.m
% Aim: use data to plot Vtu
% 2022-05-15 王良锋, 乔渭阳, 纪良, 余索远. "轴流风扇/压气机管道周向声模态的测量[J]". 航空动力学报, 2014, 29(4): 918-927.
% original coded by dgm
% modified by wjq - 2022-05-16

clc;
clear;
close all;

%% add subfunction
addpath(genpath('.'));
addpath(genpath('.'));
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

%% CPSD and phase
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
for i_file =Num_file
    eval(['load ''',chemin,'/','RotaryTest-12000-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
    Tdata{i_file}=Data(:,1:13);
    [temp_ref,freq] = cpsd(Tdata{i_file}(:,13),Tdata{i_file}(:,13),Wind,Noverlap,Nfft,Fs);
    temp_ref = sqrt(temp_ref);
    for k=1:nk
        [temp,freq] = cpsd(Tdata{i_file}(:,k),Tdata{i_file}(:,13),Wind,Noverlap,Nfft,Fs);
        CC1(:, i_file+NumSM*(k-1)) = temp./temp_ref;                                                 %360个测点的数据矩阵
    end
end

%% 获取模态信息

nk_enlarge=NumSM*nk;                                                                             %传声器可测量的模态总数
mode=-nk_enlarge/2+1 : nk_enlarge/2 -1;
for k=1:length(mode)
    a_mf(:,k)=1/nk_enlarge*CC1(:,1:nk_enlarge)*exp(mode(k)*1i*2*pi*(0:nk_enlarge-1)/nk_enlarge).';    %求解模态系数
end

%% 绘图

h=figure('Visible', 'on');
set(gcf,'position',[200 100 800 600]);
GAMMA =10*log10(abs(a_mf)/4e-10);
imagesc([-nk_enlarge/2:nk_enlarge/2],freq,GAMMA);ylim([1,rotor_speed/60*29*3.2]);
xlim([-100,100]);
axis xy;


% 切面图
Freq_slice = [1/29,1/2,1,2,3];    %对应BPF
df =freq(2) - freq(1);
h=figure('Visible', 'on');
set(gcf,'outerposition',get(0,'screensize'));%最大化
for k=1:length(Freq_slice)
    Wavemode(k,:)=max(abs(GAMMA(floor(rotor_speed/60*29*Freq_slice(k)/df)+[floor(-5/df):floor(5/df)],:)));
end
bar(mode,Wavemode');hold on
legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',12);

set(gca,'XTick',mode);
set(gca,'Ygrid','on')
    title({['模态分析'];['转速: ',num2str(rotor_speed),'-采样率：',num2str(Fs)]},'FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
ylim([80 110]);
xlim([-16 16])
