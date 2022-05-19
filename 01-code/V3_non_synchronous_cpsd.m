% Aim: 利用试验数据生成互谱矩阵
% % 输出整个频率

clc;
clear;
%close all;

%% add subfunction
addpath(genpath(['.']));
chemin = '../database/01-rotateMic';

%% add Basic parameters

zH = 0.4;         % 测试距离
NumMic = 12;          % 传声器的数量
NumSM= 30;        % 测量的次数
Radius=0.2;          % 管道半径
S=pi*Radius^2;         % 管口面积
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
    Tdata=Data(:,1:12);
    T1=  kron(ones(1,12), Tdata  );
    T2=  kron(Tdata, ones(1,12)  );
    [temp,freq]=cpsd(T1,T2,Wind,Noverlap,Nfft,Fs);
    CC1=[CC1 temp]; 
end

% 重新生成CC矩阵,注意排序顺序，reshape([1:144],24,6) 代验证
CC2=reshape(CC1,length(freq),144,30);
CC3=reshape(CC2,length(freq),12,12,30);

toc

% 切面图
Freq_slice = [1];    %对应BPF
df =freq(2) - freq(1);
f0=rotor_speed/60*29*Freq_slice;


amf=reshape(CC3(:,1,1,1),length(freq),1);
figure % 01-验证 & 用以找到最高点
plot(freq,amf);


for k=1:length(Freq_slice)
    xunhao_around=floor(f0(k)/df)+[-3:3];
    xuhao=xunhao_around(1)+find(amf(xunhao_around)==max(amf(xunhao_around)))-1;
    CC=[];
    for i_file=1:30
        CC = blkdiag(CC,reshape(CC3(xuhao,:,:,i_file),12,12));
    end
end
h=figure
imagesc(abs(CC));
axis equal


%% 传声器阵列%%%%% （柱坐标）
mic = [];
XM=Radius*ones(NumMic*NumSM,1);
YM=zeros(NumMic*NumSM,1);
for  j=1:NumSM
    theta1=[0:30:330]';
    YM((1:NumMic)+(j-1)*NumMic,1)=theta1+(j-1)*1;
end
YM=YM/180*pi;
ZM1=zH*ones(NumMic,1);
for  j=1:NumSM
    ZM((1:NumMic)+(j-1)*NumMic,1)=ZM1;
end
mic_loc=[XM,YM,ZM];



%% 构建模态系数矩阵%%%%%%%%
m = [-50:50]; %周向模态限制范围，自动删选
n = [1];      %径向模态限制范围
M =  0;       %管道流速
[G,index_mn]=matrix_G_basis(f0(1),Radius,M,mic_loc,m,n);
cond(G)

%% 非同步测量空间基函数的确定
[U,S,V] = svd(G);
Phi_basis  = U;
psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
D_measured = CC;                   % measured matrix


%% 非同步测量算法代码
%algorithm=1,2,3对应选择不同算法
algorithm=2;
switch algorithm
    case 1
       disp('FISTA');
        Li = 1.2;                       % step size of gradient descent
        N_iter = 2;                     % maximum iteation for each loop
        SC = 1e-3;                      % stopping criteria
        mu = norm(D_measured);          % starting value of continuation technique
        mu_final = 1e-16 * mu;          % final value of continuation technique
        tic
        [R_matrix_1, para_val]  = FISTA(Li,N_iter,D_measured,SC,mu_final,mu,psi_B);
        toc;
        SppSmFit = R_matrix_1;

    case 2
         disp(' ADMM');
         Omega = zeros(size(D_measured));
         Omega(find(D_measured~=0)) = 1;  % the positions which the measurements are nonzeros
         [m, n] = size(Omega);            % dimension of matrix
         SC = 0.005;                      % stopping criteria
         mIter = 14;                      % maximum iteation
         gama =2.6;       %relaxation parameter
         alpha = 28.5e-3; %regularization parameter
         mu=24.5/n;       %penalty parameter
         tic
         [R_matrix_1,err] = ADMM(D_measured, psi_B, SC, mIter, gama, mu, alpha );
          toc;
          SppSmFit = R_matrix_1;

    otherwise
         disp('CP');
         Hard_Th = 4;
         max_it = 6000;
         SC = 1e-3;
         meanORG = 0;
         varORG = 0;
         [R_matrix_1,para_AP]  = AP_cycle(Hard_Th,D_measured,psi_B,max_it,SC,meanORG,varORG);
         SppSmFit = R_matrix_1;
end

figure;
imagesc(abs(R_matrix_1))

%% 由互谱矩阵获取声压列向量
%  Spp=CC;
Spp=R_matrix_1;
P_amplitude = sqrt(diag(Spp));
P_phase = angle(Spp(:,1)/Spp(1,1));
P_complex = P_amplitude.*exp(1i*P_phase);
P=P_complex;

%% 模态识别方法:1：最小二乘法%%%%
q_re1=(G'*G)^-1*G'*P;   %03-30
abs_q1=abs(q_re1);%/(2*10-5); % pref=2e-5;
abs_q1(find(abs_q1<0))=0;
% 按照G的模态顺序重新排序
% [mode_new,order_mode]=sort(mode_prop2(:,1));
%figure; bar(mode_new,abs_q1(order_mode))
figure; bar(index_mn(:,1),abs_q1);
set(gcf,'position',[50 400 1400 300]);
set(gca,'FontSize',14)


