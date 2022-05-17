% 8000-rmp——vplot.m
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

%% CPSD and phase
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
for i_file =Num_file
    eval(['load ''',chemin,'/','RotaryTest-12000-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
    Tdata{i_file}=Data(:,1:13);
    % Step01: 通过key signal将其分段,生成多个block，每个block 6 round, 历时 25 秒。
    % 在这里需要增加等角度采样的操作：
    [key_pulse,rotor_speed]=keyRotation(Data(:,14),Fs);
    cut_number=floor(length(key_pulse)/6)-1;
    data_resample_interval(i_file)=key_pulse(7)-key_pulse((1));  
    for kb=1:cut_number
        tmp=Tdata{i_file}(key_pulse((1+(kb-1)*6)) :key_pulse(1+(kb*6)),:);
        data_block{kb,i_file}=resample(tmp,data_resample_interval(1),size(tmp,1));
    end
end
    % Step02: ensember average 得到tonal noise,历时 2.707163 秒。
      data_block_3d = reshape(cell2mat(data_block.'),data_resample_interval(1)*NumSM,13,cut_number);
      data_tonal_rms=mean(data_block_3d,3);
      data_tonal_rms2=mat2cell(data_tonal_rms,data_resample_interval(1)*ones(NumSM,1),[13]).'; % 形式与Tdata保持一致
    % Step03: r(t)=p(t)-s(t)
      data_tonal=kron(ones(cut_number,1),cell2mat(data_tonal_rms2));
      data_broadband=cell2mat(data_block)-data_tonal;

%% 作图

[Gx0,Gxx0,Fx0] = avgGxx('hann',50,'ACF',10,Fs+135,3200,Tdata{1, 1}(:,1)); %暂时fs手动微调
[Gx1,Gxx1,Fx1] = avgGxx('hann',50,'ACF',10,Fs+135,3200,data_tonal(:,1)); 
[Gx2,Gxx2,Fx2] = avgGxx('hann',50,'ACF',10,Fs+135,3200,data_broadband(:,1)); 
abs_q0=20*log10(abs(Gx0)/(2*10-5));
abs_q1=20*log10(abs(Gx1)/(2*10-5));
abs_q2=20*log10(abs(Gx2)/(2*10-5));
figure;plot(Fx0,abs_q0,'k','LineWidth',5);hold on;plot(Fx1,abs_q1,'b','LineWidth',2);   plot(Fx2,abs_q2,'r','LineWidth',2); 
xlim([0 50000])
grid on
grid minor
xlabel('Frequency/Hz')
ylabel('Sound pressure level/dB')




