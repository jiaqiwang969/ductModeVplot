% Aim: 利用试验数据生成互谱矩阵
% % 输出整个频率

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
    [temp_ref,freq] = cpsd(Data(:,13),Data(:,13),Wind,Noverlap,Nfft,Fs);
    temp_ref = sqrt(temp_ref);
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


% 切面图
Freq_slice = [1];    %对应1xBPF
df =freq(2) - freq(1);

amf=reshape(CC3(:,1,1,1),length(freq),1);
% figure % 01-验证 & 用以找到最高点
% plot(freq,amf);

f0=rotor_speed/60*29*Freq_slice;
% h=figure('Visible', 'on');
for k=1:length(Freq_slice)
    xunhao_around=floor(f0(k)/df)+[-3:3];
    xuhao=xunhao_around(1)+find(amf(xunhao_around)==max(amf(xunhao_around)))-1;
    CC=[];
    for i_file=1:30
        CC = blkdiag(CC,reshape(CC3(xuhao,:,:,i_file),12,12));
    end
end
% imagesc(abs(CC));
% axis equal


%% 选择管道内可传播的模态%%%%%%%%
% Kappa：无坐标系，r=1，有坐标系
% omega 无坐标系
f0=f0(1);
omega=2*pi*f0; % 角速度
k=omega/343; % 波数
load('Kappa.mat');    % 无流状态的管道声传播，其实是需要修正的（有流动的情况）
Kappa = Kappa/a;
Kappa=Kappa(:,1); % 只考虑周向模态

%% 传声器阵列%%%%% （柱坐标）
mic_loc = zeros(NumSM*nk,3);
XM=a*ones(nk*NumSM,1);
YM=zeros(nk*NumSM,1);
for  j=1:NumSM
    theta1=[0:30:330]';
    YM((1:nk)+(j-1)*nk,1)=theta1+(j-1)*1;
end
YM=YM/180*pi;
ZM1=0.4*ones(nk,1);
for  j=1:NumSM
    ZM((1:nk)+(j-1)*nk,1)=ZM1;
end
mic_loc=[XM,YM,ZM];


%% 根据截至频率，计算出可传播模态%
mode_prop2=propagated_models(k,Kappa);  % 可传播模态
[row,col] = size(mode_prop2);           % 可传播模态数量
[G]=matrix_G_trial(mode_prop2,Kappa,k,a,mic_loc);
%cond(G)

%% 非同步测量空间基函数的确定
[U,S,V] = svd(G);
Phi_basis  = U;
psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
D_measured = CC;                % measured matrix


%% 非同步测量算法
%  ADMM算法
Omega = zeros(size(D_measured));
Omega(find(D_measured~=0)) = 1;  % the positions which the measurements are nonzeros
[m, n] = size(Omega);            %dimension of matrix
SC = 0.005;                        % stopping criteria
mIter = 14;                       % maximum iteation
gama =2.6;      %relaxation parameter
alpha = 28.5e-3;   %regularization parameter
mu=24.5/n;       %penalty parameter
tic
[R_matrix_1,err] = ADMM(D_measured, psi_B, SC, mIter, gama, mu, alpha );
toc;

figure;
subplot(1,2,1)
imagesc(abs(CC));
axis equal
subplot(1,2,2)
imagesc(abs(R_matrix_1))
axis equal

%% 由互谱矩阵获取声压列向量
%  Spp=CC;
Spp=R_matrix_1;
P_amplitude = sqrt(diag(Spp));
P_phase = angle(Spp(:,1)/Spp(1,1));
P_complex = P_amplitude.*exp(1i*P_phase);
P=P_complex;

%% 模态识别方法:1：最小二乘法%%%%
q_re1=(G'*G)^-1*G'*P;   %03-30
pref=2e-5;
mode=mode_prop2(:,1)
abs_q1=abs(q_re1);
[I1,I2]=sort(mode)
abs_q1(I2)
figure
bar(I1,abs_q1(I2))






