clc;
clear;
% close all;

chemin = '../database/01-rotateMic';
%% 基础参数设置
Radius=0.2;        % 管道半径
rhoAir = 1.2;      % 空气密度
c = 340;           % 声速
zH = 0.4;          % 阵列的Z坐标
NumMic = 12;       % 传声器的数量
NumSM= 30;         % 非同步测量的次数
f0 = 8000/60*29*1;% 分析频率
Fs = 102400 ;      % 采样频率
time=5;            % 采样时长
Nw =  1024*5;      % length of snapshot，resolution =   Fs/Nw；
nfft =2*(Nw);      % nfft = 2^ (Nw);
F = Fs*(0:nfft-1)/nfft;   % 频率画点
w = hanning(Nw);    % window function
w = w/norm(w);      % w = w*sqrt(2/nfft); % calibration
[temp,indice_f] = min(abs(F-f0));       %找到最接近波束形成频率f的频率所在列数
Nsnap = 100;        % Nsnap = floor(2*((Fs*10))/Nw-1);%快照数量


%% 构建互谱矩阵
pM_temp = zeros(NumMic,Nsnap);
pM = [];
p=zeros(Fs*time,NumMic);
Ind = [1:NumSM];
Num_file = Ind ;

for k =Num_file
    eval(['load ''',chemin,'/','RotaryTest-8000-Rotate-No-',num2str(k),'.mat''']);
    for index_M = 1:NumMic
        p(:,index_M)=Data(:,index_M);
        for m = 1:Nsnap
            Pf = conj(fft(p((m-1)*Nw/2+(1:Nw),index_M).*w,nfft));
            pM_temp(index_M,m) = Pf(indice_f);
        end
    end
    pM = [pM ; pM_temp];
end

% Cross spectral matrix computation (CSM matrix): data missing matrix
MM = zeros(NumMic*NumSM,NumMic*NumSM);
Cov_pM = zeros(NumMic,NumMic,NumSM);
for n = 1:NumSM
    Cov_pM(:,:,n) = pM((n-1)*NumMic+1:n*NumMic,:)*pM((n-1)*NumMic+1:n*NumMic,:)'/Nsnap;
    Cov_pM(:,:,n) = (squeeze(Cov_pM(:,:,n)) + squeeze(Cov_pM(:,:,n))')/2;%squeeze函数用于删除矩阵中的单一维
end

for n = 1:NumSM
    MM(((n-1)*NumMic+(1:NumMic)),((n-1)*NumMic+(1:NumMic))) = squeeze(Cov_pM(:,:,n));
end

figure
imagesc(abs(MM));



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
[G,index_mn]=matrix_G_basis(f0,Radius,M,mic_loc,m,n);
cond(G)

% omega=2*pi*f0; % 角速度
% k=omega/343;   % 波数
% load('Kappa.mat');    % 无流状态的管道声传播，其实是需要修正的（有流动的情况）
% Kappa = Kappa/Radius;
% Kappa=Kappa(:,1);     % 只考虑周向模态
% mode_prop2=propagated_models(k,Kappa);  % 可传播模态
% [mode_new,order_mode]=sort(mode_prop2(:,1))
% [G]=matrix_G_trial(mode_prop2,Kappa,k,Radius,mic_loc);
% cond(G)


%% 非同步测量空间基函数的确定
[U,S,V] = svd(G);
Phi_basis  = U;
psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
D_measured = MM;                   % measured matrix


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
%figure; bar(mode_new,abs_q1(order_mode))
figure; bar(index_mn(:,1),abs_q1);
set(gcf,'position',[50 400 1400 300]);
set(gca,'FontSize',14)


