% Aim: 利用试验数据生成互谱矩阵
% % 输出整个频率

clc;
clear;
close all;

%% add subfunction
addpath(genpath(['.']));
chemin = '../database/01-rotateMic';

%% add Basic parameters

zH = 0.4;           % 测试距离
zRef=zH-0.08;      % 参考传声器到测试面的轴向距离
NumMic = 12;        % 传声器的数量
NumSM  = 30;        % 测量的次数
Radius = 0.2;       % 管道半径
Area = pi*Radius^2; % 管口面积
Fs = 102400;        % 采样频率
time = 5;           % 采样时间

rotor_speed=14000;  % 轴转速信息

%% data processing
L_signal = Fs*time;             %信号长度
L_seg = round(L_signal/100);    %确定对信号处理划窗长度
Wind = hamming(L_seg);          %确定对数据进行汉宁窗处理
Noverlap = round(L_seg/2);      %确定信号划窗重叠率
Nfft = 2^(ceil(log2(L_seg))+1); %确定分析频率
round=7;                        %分段,生成多个block，每个block 为round

%% data processing
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
for i_file =Num_file
    eval(['load ''',chemin,'/','RotaryTest-',num2str(rotor_speed),'-Rotate-No-',num2str(i_file),'.mat''']);       %读取数据
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

%% 是否将参考传声器植入到互功率谱矩阵中
if_consider_ref=1 %with method1: reference sensor is same to others
%% CPSD and phase
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
CC1=[];
for i_file = Num_file
    temp_data=data_all(:,(1+(i_file-1)*13):(1+(i_file)*13-2+if_consider_ref)); % 把传声器也考虑在内！！
    T1=  kron(ones(1,12+if_consider_ref), temp_data  );
    T2=  kron(temp_data, ones(1,12+if_consider_ref)  );
    [temp,freq]=cpsd(T1,T2,Wind,Noverlap,Nfft,Fs);
    CC1=[CC1 temp];
end
% 重新生成CC矩阵,注意排序顺序，reshape([1:144],24,6) 代验证
CC2=reshape(CC1,length(freq),(12+if_consider_ref)*(12+if_consider_ref),NumSM);
CC3=reshape(CC2,length(freq),12+if_consider_ref,12+if_consider_ref,NumSM);
amf=reshape(CC3(:,1,1,1),length(freq),1);
% figure % 01-验证 & 用以找到最高点
% plot(freq,amf);
%
Freq_slice = [1 2];    %对应1xBPF
df =freq(2) - freq(1);
f0=rotor_speed/60*29*Freq_slice;


for k=1:length(Freq_slice)
    xunhao_around=floor(f0(k)/df)+[-3:3];
    xuhao=xunhao_around(1)+find(amf(xunhao_around)==max(amf(xunhao_around)))-1;
    CC{k}=[];
    for i_file=1:NumSM
        CC{k} = blkdiag(CC{k},reshape(CC3(xuhao,:,:,i_file),12+if_consider_ref,12+if_consider_ref));
    end
    % 对CC进行处理，形成真正对ref参考矩阵
    if if_consider_ref==1
        % Step01: 裁掉传声器所在对行
        CC_1=CC{k};
        CC_1(13*[1:30],:)=[];
        CC_1(:,13*[1:30])=[];
        % Step02: 补上参考传声器的影响
        CC_2=CC{k};CC_2(13*[1:30],:)=[]; temp1=sum(CC_2(:,13*[1:30]),2);
        CC_3=CC{k};CC_3(:,13*[1:30])=[]; temp2=sum(CC_3(13*[1:30],:),1);
        CC_diag=diag(CC{k});temp3=mean(CC_diag(13*[1:30]));
        % Step03: 重新组装
        CC{k}=[[temp3 temp2];[temp1 CC_1]];
    end
end 




figure;
subplot(1,2,1)
imagesc(abs(CC{1}));axis equal %共轭对称矩阵
subplot(1,2,2)
imagesc(abs(CC{2}));axis equal %共轭对称矩阵




%% 传声器阵列%%%%% （柱坐标）
mic_loc = zeros(NumSM*NumMic,3);
XM=Radius*ones(NumSM*NumMic,1);
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

if if_consider_ref==1
    mic_loc=[[Radius 0 zRef];mic_loc];
end

figure
xs = mic_loc(:,1).*cos(mic_loc(:,2));ys=mic_loc(:,1).*sin(mic_loc(:,2));
for i = 1:NumSM
    plot3(xs((1:NumMic)+(i-1)*NumMic),ys((1:NumMic)+(i-1)*NumMic),mic_loc((1:NumMic)+(i-1)*NumMic,3),'.');
    hold on;
end

%% 根据截至频率，计算出可传播模态
% omega=2*pi*f0; % 角速度
% k=omega/343;   % 波数
% load('Kappa.mat');    % 无流状态的管道声传播，其实是需要修正的（有流动的情况）
% Kappa = Kappa/Radius;
% Kappa=Kappa(:,1);     % 只考虑周向模态
% mode_prop2=propagated_models(k,Kappa);  % 可传播模态
% [G1]=matrix_G_trial(mode_prop2,Kappa,k,Radius,mic_loc);
% cond(G1)



for k=1:length(Freq_slice)

    %% 该处替换为wjq写的green函数
    % 输入：w,Radius,mic_loc
    % 输出：G:361*[-mn,+mn] 频率域

    m = [-50:50];   % 周向模态限制范围，自动删选
    n = [1];        % 径向模态限制范围
    M =  0.1;       % 管道流速
    [G{k},index_mn{k}]=matrix_G_basis(f0(k),Radius,M,mic_loc,m,n);
    cond(G{k})


    %% 非同步测量空间基函数的确定
    [U,S,V] = svd(G{k});
    Phi_basis  = U;
    psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
    D_measured = CC{k};                % measured matrix

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


    % figure; subplot(1,2,1);imagesc(abs(CC));axis equal
    % subplot(1,2,2);imagesc(abs(R_matrix_1));axis equal

    %% 由互谱矩阵获取声压列向量
    %  Spp=CC;
    Spp=R_matrix_1;
    P_amplitude = sqrt(diag(Spp));
    P_phase = angle(Spp(:,1)/Spp(1,1));
    P_complex = P_amplitude.*exp(1i*P_phase);
    P=P_complex;

    %% 模态识别方法:1：最小二乘法%%%%
    q_re1{k}=(G{k}'*G{k})^-1*G{k}'*P;   %03-30
end


figure; 
bar(index_mn{1}(:,1),abs(q_re1{1}));
hold on;
bar(index_mn{2}(:,1),abs(q_re1{2}));
legend({'1*BPF';'2*BPF';},'Location','NorthEast','FontSize',12);








