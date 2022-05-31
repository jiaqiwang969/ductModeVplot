% Aim: 利用试验数据生成互谱矩阵
% % 输出整个频率

clc;
clear;
close all;

%% add subfunction
addpath(genpath(['../']));
chemin = '../../database/01-rotateMic';

%% add Basic parameters

zH = 0;             % 测试距离
zRef=zH-0.08;       % 参考传声器到测试面的轴向距离
NumMic = 12;        % 传声器的数量
NumSM  = 30;        % 测量的次数
Radius = 0.2;       % 管道半径
Area = pi*Radius^2; % 管口面积
Fs = 102400;        % 采样频率
time = 5;           % 采样时间

rotor_speed=14000;   % 轴转速信息

%w=2*pi*rotor_speed/60*29/343*Radius;

%% data processing
L_signal = Fs*time;             %信号长度
L_seg = round(L_signal/100);    %确定对信号处理划窗长度
Wind =  hamming(L_seg);         %确定对数据进行汉宁窗处理
Noverlap = round(L_seg/2);      %确定信号划窗重叠率
Nfft = 2^(ceil(log2(L_seg))+1); %确定分析频率
block=2;%分段,生成多个block，每个block 为round



%% data processing-resample
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
for i_file =Num_file
    eval(['load ''',chemin,'/','RotaryTest-',num2str(rotor_speed),'-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
    Tdata{i_file}=Data(1:L_signal,1:13);
    % Step01: 通过key signal将其分段,生成多个block，每个block 为 round, 历时 25 秒。
    % 在这里需要增加等角度采样的操作：
    [key_pulse,rotor_speed]=keyRotation(Data(1:L_signal,14),Fs);
    cut_number(i_file)=floor(length(key_pulse)/block)-1;
    data_resample_interval(i_file)=key_pulse(block+1)-key_pulse((1));
    for kb=1:cut_number(1)
        tmp=Tdata{i_file}(key_pulse((1+(kb-1)*block)):key_pulse(1+(kb*block)),:);
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
y_tonal=kron(ones(cut_number,1),cell2mat(data_tonal_rms2));
y_residual=data_all-y_tonal;





% %%  循环维纳滤波（提取二阶循环平稳信号）
% L1 = size(y_residual,1);
% Na = 10; % 循环频率的阶数
% orders = (1:Na/2); % 循环频率的个数
% orders = [orders;-orders];
% orders = orders(:); % 总的循环频率个数
% alpha = 1230/Fs; % 第一阶循环频率 /136.8
% phi = 2*pi*alpha*(0:L1-1)'; % 相位的索引
% fi = alpha*ones(L1,1); % 频率的索引
% 
% 
% tic
% parfor k=1:size(y_residual,2)
%     k
%     % 估计循环互谱矩阵
%     Nw = 2^8;
%     [Syya,Syaya] = CyclicSpecMat(y_residual(:,k),phi,fi,orders,Nw);
%     f = Syaya.f*Fs;
%     % % 构建循环维纳滤波器
%     Ns = 1; % 循环平稳声源的个数
%     [G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,1e-6);
%     % 基于Gabor变换进行滤波
%     [xhat{k},cx{k}] = CyclicSpatialFiltering(y_residual(:,k),G,phi,fi,orders);
% end
% toc
% y_residual0=cell2mat(xhat);

%% 作图
[Gx0,Gxx0,Fx0] = avgGxx('hann',50,'ACF',10,Fs,3200,data_all(:,1)); 
[Gx1,Gxx1,Fx1] = avgGxx('hann',50,'ACF',10,Fs,3200,y_tonal(:,1)); 
[Gx2,Gxx2,Fx2] = avgGxx('hann',50,'ACF',10,Fs,3200,y_residual(:,1)); 
% [Gx3,Gxx3,Fx3] = avgGxx('hann',50,'ACF',10,Fs,3200,y_residual0(:,1)); 

abs_q0=20*log10(abs(Gx0)/(20*1e-6));
abs_q1=20*log10(abs(Gx1)/(20*1e-6));
abs_q2=20*log10(abs(Gx2)/(20*1e-6));
% abs_q3=20*log10(abs(Gx3)/(20*1e-6));

figure;semilogx(Fx0,abs_q0,'k','LineWidth',5);hold on;plot(Fx1,abs_q1,'b','LineWidth',2);   plot(Fx2,abs_q2,'r','LineWidth',2); %plot(Fx3,abs_q3,'r','LineWidth',2); 
xlim([1000 100000])
ylim([50 200])
grid on
grid minor
xlabel('Frequency/Hz','FontSize',20)
ylabel('Sound pressure level/dB', 'FontSize',20)

hold on
fr=[1:29]*(14000/60);
plot([fr;fr],[50*ones(1,length(fr));200*ones(1,length(fr))],'--')
legend({'total';'CS1,单音噪声';'CS2,宽带噪声';},'Location','NorthEast','FontSize',12);



%% 方法01:CC算法
[GAMMA_all,freq_all]=CC(y_residual);
%[GAMMA_tonal,freq_tonal]=CC(y_tonal);


%% 方法02:合成孔径算法
algorithm=3;
Mode1=[-60:60];
Freq_slice_all = 0.1:0.0125/2:2.2;    
[Gamma_ref_all]=SA(y_residual,Freq_slice_all,rotor_speed,algorithm);
%Freq_slice_tonal= 1;    
%[Gamma_ref_tonal]=SA(y_tonal,Freq_slice_tonal,rotor_speed,algorithm);


%%
h=figure;
subplot(2,2,1)
bar(Mode1,GAMMA_all(find(freq_all<(rotor_speed/60*29+3)&freq_all>(rotor_speed/60*29-3)),:))
axis xy
grid minor
% colorbar
xlim([-31,31])
title('Tonal')
subplot(2,2,2)
imagesc(Mode1,freq_all/(rotor_speed/60*29),log(GAMMA_all))
ylim([0.1,2.2])
axis xy
colorbar
xlim([-50,50])
title('All Band')
subplot(2,2,3)
bar(Mode1,Gamma_ref_all(find(Freq_slice_all==1),:))
axis xy
grid minor
% colorbar
xlim([-31,31])
subplot(2,2,4)
imagesc(Mode1,Freq_slice_all,log(Gamma_ref_all))
axis xy
colorbar
ylim([0.1,2.2])
xlim([-50,50])


saveas(h,[date,'-',num2str(rotor_speed),'-',num2str(algorithm),'.fig'])
saveas(h,[date,'-',num2str(rotor_speed),'-',num2str(algorithm),'.eps'])




