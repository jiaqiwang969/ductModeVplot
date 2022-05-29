% Aim: 利用试验数据生成互谱矩阵
% % 输出整个频率

% clc;
% clear;
% close all;

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

rotor_speed=12000;              %轴转速信息
%% 为了提高计算效率，选择先降采样
%Fs_new=rotor_speed/60*29*2*2.56;

%% data processing
L_signal = Fs*time;             %信号长度
L_seg = round(L_signal/100);    %确定对信号处理划窗长度
Wind = hamming(L_seg);          %确定对数据进行汉宁窗处理
Noverlap = round(L_seg/2);      %确定信号划窗重叠率
Nfft = 2^(ceil(log2(L_seg))+1); %确定分析频率


tic
%% CPSD and phase
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
CC1=[];
for i_file =Num_file
    eval(['load ''',chemin,'/','RotaryTest-12000-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
    %Tdata=resample(Data(:,1:13),Fs,Fs_new);
    Tdata=Data(:,1:12);
    %[temp_ref,freq] = cpsd(Data(:,13),Data(:,13),Wind,Noverlap,Nfft,Fs);
    %temp_ref = sqrt(temp_ref);
    % CPSD 矩阵形式，加速运算 历时 6.863685 秒，for-loop：历时 27.389531 秒。
    T1=  kron(ones(1,12), Tdata  );
    T2=  kron(Tdata, ones(1,12)  );
    [temp,freq]=cpsd(T1,T2,Wind,Noverlap,Nfft,Fs);
    
    % figure;
    % plot(freq,abs(temp(:,5)))
    CC1=[CC1 temp];  %"./temp_ref" for tonal noise or not
end

% 重新生成CC矩阵,注意排序顺序，reshape([1:144],24,6) 代验证
CC2=reshape(CC1,length(freq),144,30);
CC3=reshape(CC2,length(freq),12,12,30);

toc

% 切面图
Freq_slice = [1];    %对应BPF
df =freq(2) - freq(1);

amf=reshape(CC3(:,1,1,1),length(freq),1);
figure % 01-验证 & 用以找到最高点
plot(freq,amf);


h=figure('Visible', 'on');
for k=1:length(Freq_slice)
    xunhao_around=floor(rotor_speed/60*29*Freq_slice(k)/df)+[-3:3];
    xuhao=xunhao_around(1)+find(amf(xunhao_around)==max(amf(xunhao_around)))-1;
    CC=[];
    for i_file=1:30
        CC = blkdiag(CC,reshape(CC3(xuhao,:,:,i_file),12,12));
    end
end
imagesc(abs(CC));axis equal;colorbar






