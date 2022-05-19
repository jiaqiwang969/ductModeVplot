% Aim: use data to plot Vtu
% 2022-05-15 
% Ref: Behn M, Pardowitz B, Tapken U. Separation of tonal and broadband 
% noise components by cyclostationary analysis of the modal sound field in 
% a low-speed fan test rig[C]//International Conference of Fan Noise, 
% Aerodynamics, Applications and Systems. 2018: 18-20.
% wjq - 2022-05-17
% - duct mode Vplot with resample
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

rotor_speed=10000;              %轴转速信息
round=7;                        %分段,生成多个block，每个block 为round

%% data processing
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
for i_file =Num_file
    eval(['load ''',chemin,'/','RotaryTest-10000-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
    Tdata{i_file}=Data(:,1:13);
    % Step01: 通过key signal将其分段,生成多个block，每个block 为 round, 历时 25 秒。
    % 在这里需要增加等角度采样的操作：
    [key_pulse,rotor_speed]=keyRotation(Data(:,14),Fs);
    cut_number(i_file)=floor(length(key_pulse)/round)-1;
    data_resample_interval(i_file)=key_pulse(round+1)-key_pulse((1));  
    for kb=1:cut_number(1)
        tmp=Tdata{i_file}(key_pulse((1+(kb-1)*round)):key_pulse(1+(kb*round)),:);
        data_block{kb,i_file}=resample(tmp,data_resample_interval(1),size(tmp,1));
    end
end
      cut_number=cut_number(1);
    % Step01-1: Data_all with resample
      data_all= cell2mat(data_block);
    % Step02: ensember average 得到tonal noise, 历时 2.707163 秒。
      data_block_3d = reshape(cell2mat(data_block.'),data_resample_interval(1)*NumSM,13,cut_number);
      data_tonal_rms=mean(data_block_3d,3);
      data_tonal_rms2=mat2cell(data_tonal_rms,data_resample_interval(1)*ones(NumSM,1),[13]).'; % 形式与Tdata保持一致
    % Step03: r(t)=p(t)-s(t)
      data_tonal=kron(ones(cut_number,1),cell2mat(data_tonal_rms2));
      data_broadband=cell2mat(data_block)-data_tonal;
    



%% 作图
[Gx0,Gxx0,Fx0] = avgGxx('hann',50,'ACF',10,Fs,3200,data_all(:,1)); 
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




%% CPSD and phase
for i_file =1:30
    temp_data_all=data_all(:,(1+(i_file-1)*13):(1+(i_file)*13-1));
    temp_data_tonal=data_tonal(:,(1+(i_file-1)*13):(1+(i_file)*13-1));
    temp_data_broadband=data_broadband(:,(1+(i_file-1)*13):(1+(i_file)*13-1));

    for k=1:nk
        [temp_all,freq] = cpsd(temp_data_all(:,k),temp_data_all(:,13),Wind,Noverlap,Nfft,Fs);
        [temp_tonal,freq] = cpsd(temp_data_tonal(:,k),temp_data_tonal(:,13),Wind,Noverlap,Nfft,Fs);
        [temp_broadband,freq] = cpsd(temp_data_broadband(:,k),temp_data_broadband(:,13),Wind,Noverlap,Nfft,Fs);
        CC1_all(:, i_file+NumSM*(k-1)) = temp_all;   
        CC1_tonal(:, i_file+NumSM*(k-1)) = temp_tonal;   
        CC1_broadband(:, i_file+NumSM*(k-1)) = temp_broadband;  
    end
end


%% 获取模态信息

nk_enlarge=NumSM*nk;                                                                             %传声器可测量的模态总数
mode=-nk_enlarge/2+1 : nk_enlarge/2 -1;
for k=1:length(mode)
    a_mf_all(:,k)=1/nk_enlarge*CC1_all(:,1:nk_enlarge)*exp(mode(k)*1i*2*pi*(0:nk_enlarge-1)/nk_enlarge).';    %求解模态系数
    a_mf_tonal(:,k)=1/nk_enlarge*CC1_tonal(:,1:nk_enlarge)*exp(mode(k)*1i*2*pi*(0:nk_enlarge-1)/nk_enlarge).';    %求解模态系数
    a_mf_broadband(:,k)=1/nk_enlarge*CC1_broadband(:,1:nk_enlarge)*exp(mode(k)*1i*2*pi*(0:nk_enlarge-1)/nk_enlarge).';    %求解模态系数

end

%% 绘图

h=figure('Visible', 'on');
set(gcf,'position',[200 100 800 600]);
subplot(2,3,1)
GAMMA_all =20*log10(abs(a_mf_all)/(2*10-5));
imagesc([-nk_enlarge/2:nk_enlarge/2],freq,GAMMA_all);ylim([1,rotor_speed/60*29*3.2]);
xlim([-100,100]);
axis xy;
title("all noise")
subplot(2,3,2)
GAMMA_tonal =20*log10(abs(a_mf_tonal)/(2*10-5));
imagesc([-nk_enlarge/2:nk_enlarge/2],freq,GAMMA_tonal);ylim([1,rotor_speed/60*29*3.2]);
xlim([-100,100]);
axis xy;
title("tonal noise")
subplot(2,3,3)
GAMMA_broadband =20*log10(abs(a_mf_broadband)/(2*10-5));
imagesc([-nk_enlarge/2:nk_enlarge/2],freq,GAMMA_broadband);ylim([1,rotor_speed/60*29*3.2]);
xlim([-100,100]);
axis xy;
title("broadband noise")



% 切面图
Freq_slice = [1,2];    %对应BPF
df =freq(2) - freq(1);
subplot(2,3,4)
for k=1:length(Freq_slice)
    Wavemode_all(k,:)=max(abs(a_mf_all(floor(rotor_speed/60*29*Freq_slice(k)/df)+[floor(-5/df):floor(5/df)],:)));
end
bar(mode,Wavemode_all');hold on
legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',12);
set(gca,'XTick',mode);
set(gca,'Ygrid','on')
    title({['模态分析'];['转速: ',num2str(rotor_speed),'-采样率：',num2str(Fs)]},'FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
% ylim([80 110]);
xlim([-16 16])
subplot(2,3,5)
for k=1:length(Freq_slice)
    Wavemode_tonal(k,:)=max(abs(a_mf_tonal(floor(rotor_speed/60*29*Freq_slice(k)/df)+[floor(-5/df):floor(5/df)],:)));
end
bar(mode,Wavemode_tonal');hold on
legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',12);
set(gca,'XTick',mode);
set(gca,'Ygrid','on')
    title({['模态分析'];['转速: ',num2str(rotor_speed),'-采样率：',num2str(Fs)]},'FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
% ylim([80 110]);
xlim([-16 16])
subplot(2,3,6)
for k=1:length(Freq_slice)
    Wavemode_broadband(k,:)=max(abs(a_mf_broadband(floor(rotor_speed/60*29*Freq_slice(k)/df)+[floor(-5/df):floor(5/df)],:)));
end
bar(mode,Wavemode_broadband');hold on
legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',12);
set(gca,'XTick',mode);
set(gca,'Ygrid','on')
    title({['模态分析'];['转速: ',num2str(rotor_speed),'-采样率：',num2str(Fs)]},'FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
% ylim([80 110]);
xlim([-16 16])



