% 8000-rmp——vplot.m
% Aim: use data to plot Vtu
% 2022-05-11 ymj

clc;
clear;
close all;

%% add subfunction
addpath(genpath('E:\模态识别\基于非同步测量的压气机管道模态识别数据\Testcode-ver1\'));
addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\simulation\'));
addpath(genpath('E:\模态识别\模态识别代码\mode_detect\mode_detect\simulation\subprogram\'));
chemin = 'E:\模态识别\基于非同步测量的压气机管道模态识别数据\实验20-2019-11-16旋转机闸 测试';

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
    
    rotor_speed=8000;               %轴转速信息
    
    %% CPSD and phase 
        Ind = [1:NumSM];   %设定循环次数
         Num_file = Ind ;         
       for i_file =Num_file
    eval(['load ''',chemin,'\','RotaryTest-8000-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
       Tdata{i_file}=Data(:,1:13);
       phase_info{i_file}=cpsd(Tdata{i_file}(:,13),Tdata{1}(:,13),Wind,Noverlap,Nfft,Fs)./...
           abs(cpsd(Tdata{i_file}(:,13),Tdata{1}(:,13),Wind,Noverlap,Nfft,Fs));                 %获取相位
          for k=1:nk
            [temp,freq] = cpsd(Tdata{i_file}(:,k),Tdata{i_file}(:,13),Wind,Noverlap,Nfft,Fs);   %CPSD
              CC1(:,i_file+NumSM*(k-1)) = temp;                                                 %360个测点的数据矩阵
                                 %            ./phase_info{i_file};                             %相位修正
        end
 end

 %% 获取模态信息

    nk_enlarge=NumSM*nk;                                                                             %传声器可测量的模态总数
    m=-nk_enlarge/2:nk_enlarge/2;            
    for k=1:length(m)
        a_mf(:,k)=1/nk_enlarge*CC1(:,1:nk_enlarge)*exp(m(k)*1i*2*pi*(1:nk_enlarge)/nk_enlarge).';    %求解模态系数
    end

    %% 绘图

      h=figure('Visible', 'on');
      set(gcf,'position',[200 100 800 600]);
     GAMMA =10*log10(abs(a_mf)/4e-10);
     imagesc([-nk_enlarge/2:nk_enlarge/2],freq,GAMMA);ylim([1,rotor_speed/60*29*3.2]);
       xlim([-100,100]); 
     axis xy;