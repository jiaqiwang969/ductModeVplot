clear;clc;
close all
addpath(genpath(['.']));
%% initial parameters 初始参数%
fa=4833; % 分析频率x
omega=2*pi*fa; % 角速度
r=0.2;         % 半径
c=340; k=omega/c; rho=1.225;pref=4e-10;
S=pi*r^2;        % 面积
zH = 0;       % Z coordinate of the array (of the measurements)
M = 12;     % number of microphones
Nm= 30;        % number of sequential measurements
Nsnap = 100;       % number of snapshots

%% 传声器阵列%%%%% （柱坐标）
        mic = [];
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
        mic=[XM,YM,ZM];
        


%% The propagation mode numbers 根据截至频率，计算出可传播模态%
load('Kappa.mat');
Kappa = Kappa/r; 
Kappa=Kappa(:,1); % 只考虑周向模态
mode_prop2=propagated_models(k,Kappa);  % 可传播模态
[row,col] = size(mode_prop2); % 可传播模态数量        
    %确定目标模态
      m_p = [3,6]; n_p = [0];    %%需要产生的模态阶数
      AA=zeros(1,row);
      for ii=m_p
      ind = ismember(mode_prop2,[ii,n_p],'rows')';
      posi = find(ind==1);
      AA(posi)=1;                  %%%%%目标模态系数为1，其他模态系数为0%%
      end
%% 计算麦克风处的声压值%
       [P1]=presure1(mic,mode_prop2,Kappa,k,r,AA);
       P1 =P1';
        a = normrnd(1,1,1,Nsnap);
        Pdirect=P1*a;
        
        %信噪比控制
        SNR=30;
        pM = AddNoise(Pdirect,SNR,Nm,M,Nsnap);     %         对声压信号添加噪声
             
          Spp1 = pM*pM'/Nsnap;
          SppOrg = (Spp1 + Spp1')./2;             % 相当于用 NumSM*Single 个麦克风进行测量而获得的矩阵(加噪声) 
          figure;
          imagesc(abs(SppOrg));colormap(flipud(hot));
          colorbar;
          axis square;xlabel('Microphone Index');ylabel('Microphone Index');
          set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);       



        %% csm构建
        % Matrice de covariance des mesures
Cov_pM = zeros(M,M,Nm);
for m = 1:Nm
    Cov_pM(:,:,m) = pM((m-1)*M+(1:M),:)*pM((m-1)*M+(1:M),:)'/Nsnap;
end

MM = zeros(Nm*M);

for m = 1:Nm
    MM(((m-1)*M+(1:M)),((m-1)*M+(1:M))) = squeeze(Cov_pM(:,:,m));
end
figure,imagesc(abs(MM))
        
%% 构建观测矩阵%%%%%%%%
         [G]=matrix_G_trial(mode_prop2,Kappa,k,r,mic);
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
%         mu_final = 1e-16 * mu;           % final value of continuation technique
%         tic
%         [R_matrix_1, para_val]  = FISTA(Li,N_iter,D_measured,SC,mu_final,mu,psi_B);
%         toc;
%         SppSmFit = R_matrix_1;      
% 
%         
        %%   ADMM

 
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
P_complex = P_amplitude.*exp(1i*P_phase);
P=P_complex;

 figure;
imagesc(abs(Spp));colormap(flipud(hot));%caxis([0,2500]);
colorbar;
axis square;xlabel('Microphone Index');ylabel('Microphone Index');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);  


%% 方法:1：最小二乘法%%%%
q_re1=(G'*G)^-1*G'*P;   %03-30
abs_q1=(abs(q_re1));%/(2*10-5); % pref=2e-5;
abs_q1(find(abs_q1<0))=0;
% 按照G的模态顺序重新排序
[mode_new,order_mode]=sort(mode_prop2(:,1));
figure; bar(mode_new,abs_q1(order_mode))
set(gcf,'position',[50 400 1400 300]);
set(gca,'FontSize',14)