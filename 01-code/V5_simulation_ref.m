clear;clc;
close all
addpath(genpath(['.']));
%% initial parameters 初始参数%
fa=4833; % 分析频率
omega=2*pi*fa; % 角速度
r=0.2;         % 半径
c=343; k=omega/c; rho=1.225;pref=2e-5;

S=pi*r^2;      % 面积
zH = 0.4;      % Z coordinate of the array (of the measurements)
zH_ref=zH-0.08;
M = 12;        % number of microphones
Nm= 30;        % number of sequential measurements
Nsnap = 100;   % number of snapshots
N_ref=1;       % 参考传声器的数量

%% 传声器阵列%%%%% （柱坐标）
mic_loc = [];
XM=r*ones(M*Nm,1);
YM=zeros(M*Nm,1);
for  j=1:Nm
    theta1=[0:30:330]';
    YM((1:M)+(j-1)*M,1)=theta1+(j-1)*1;
end
YM=YM/180*pi;
ZM1=zH*ones(M,1);
for  j=1:Nm
    ZM((1:M)+(j-1)*M,1)=ZM1;
end
mic_loc=[XM,YM,ZM];
ref_mic=[r,0,zH_ref];
mic=[ref_mic;mic_loc];


%% The propagation mode numbers 根据截至频率，计算出可传播模态%
load('Kappa.mat');
Kappa = Kappa/r;
Kappa=Kappa(:,1); % 只考虑周向模态
mode_prop2=propagated_models(k,Kappa);  % 可传播模态
[row,col] = size(mode_prop2); % 可传播模态数量
%确定目标模态
m_p = [-3,6]; n_p = [0];    %%需要产生的模态阶数
AA=zeros(1,row);
for ii=m_p
    ind = ismember(mode_prop2,[ii,n_p],'rows')';
    posi = find(ind==1);
    AA(posi)=1;                  %%%%%目标模态系数为1，其他模态系数为0%%
end

%% 计算麦克风处的声压值%
[P1]=presure1(mic_loc,mode_prop2,Kappa,k,r,AA);
P1 =P1';
P_ref=presure1(ref_mic,mode_prop2,Kappa,k,r,AA);
a = normrnd(1,1,1,Nsnap);
Pdirect=P1*a;
P_ref=P_ref*a;
%信噪比控制
SNR=30;
pM = AddNoise(Pdirect,SNR,Nm,M,Nsnap);  % 对声压信号添加噪声
P_ref=AddNoise(P_ref,SNR,1,1,Nsnap);

Spp1 = pM*pM'/Nsnap;
SppOrg = (Spp1 + Spp1')./2;             % 相当于用 NumSM*Single 个麦克风进行测量而获得的矩阵(加噪声)
figure;
imagesc(abs(SppOrg));colormap(flipud(hot));
colorbar;
axis square;xlabel('Microphone Index');ylabel('Microphone Index');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);


%% csm构建
Cov_pM = zeros(M,M,Nm);
for m = 1:Nm
    Cov_pM(:,:,m) = pM((m-1)*M+(1:M),:)*pM((m-1)*M+(1:M),:)'/Nsnap;
end

%  references
Cov_ref = P_ref*P_ref'/Nsnap;

% 边界参照矩阵生成
Cov_ref_pM = zeros(N_ref,M,Nm);
for m = 1:Nm
    Cov_ref_pM(:,:,m) = P_ref*pM((m-1)*M+(1:M),:)'/Nsnap;
end

MM = zeros(Nm*M+N_ref);
MM(1:N_ref,1:N_ref) = Cov_ref;
for m = 1:Nm
    MM(1:N_ref,N_ref+(m-1)*M+(1:M)) = squeeze(Cov_ref_pM(:,:,m));
    MM(N_ref+((m-1)*M+(1:M)),1:N_ref) = squeeze(Cov_ref_pM(:,:,m))';
end
for m = 1:Nm
    MM(N_ref+((m-1)*M+(1:M)),N_ref+((m-1)*M+(1:M))) = squeeze(Cov_pM(:,:,m));
end

% Original matrix
Cov_all = pM*pM'/Nsnap;
ORGmatrixdataX = [MM(1:N_ref,:);[MM(N_ref+1:Nm*M+N_ref,1:N_ref),Cov_all]];

figure,imagesc(abs(ORGmatrixdataX))
figure,imagesc(abs(MM))

%% 构建观测矩阵%%%%%%%%
%[G]=matrix_G_trial(mode_prop2,Kappa,k,r,mic);
m = [-50:50]; %周向模态限制范围，自动删选
n = [1];      %径向模态限制范围
M =  0;       %管道流速
[G,index_mn]=matrix_G_basis(fa,r,M,mic,m,n);
cond(G)
%%  获取非同步测量空间基函数
[U,S,V] = svd(G);
Phi_basis  = U;
psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
D_measured = MM;          % measured matrix
%%   fista
%         Li = 1.2;                       % step size of gradient descent
%         N_iter = 2;                     % maximum iteation for each loop
%         SC = 1e-3;                      % stopping criteria
%         mu = norm(D_measured);          % starting value of continuation technique
%         mu_final = 1e-16 * mu;          % final value of continuation technique
%         tic
%         [R_matrix_1, para_val]  = FISTA(Li,N_iter,D_measured,SC,mu_final,mu,psi_B);
%         toc;
%         SppSmFit = R_matrix_1;
%
%
%%   ADMM


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


%%  CP
% Hard_Th = 4;
% max_it = 6000;
% SC = 1e-3;
% meanORG = 0;
% varORG = 0;
% [R_matrix_1,para_AP]  = AP_cycle(Hard_Th,D_measured,psi_B,max_it,SC,meanORG,varORG);
%         SppSmFit = R_matrix_1;

%% 模态识别声压列向量
Spp=SppSmFit;
P_amplitude = sqrt(diag(Spp)); % 获取幅值
P_phase = angle(Spp(:,1)/Spp(1,1));  % 获取相位角
P = P_amplitude.*exp(1i*P_phase);

figure;
imagesc(abs(Spp));colormap(flipud(hot));
colorbar;
axis square;xlabel('Microphone Index');ylabel('Microphone Index');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);


%% 方法:1：最小二乘法%%%%
q_re1=(G'*G)^-1*G'*P;   %03-30
abs_q1=abs(q_re1);%/(2*10-5); % pref=2e-5;
abs_q1(find(abs_q1<0))=0;
% 按照G的模态顺序重新排序
% [mode_new,order_mode]=sort(mode_prop2(:,1));
%figure; bar(mode_new,abs_q1(order_mode))
figure; bar(index_mn(:,1),abs_q1);
set(gcf,'position',[50 400 1400 300]);
set(gca,'FontSize',14)