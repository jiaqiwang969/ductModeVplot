clear;clc;
close all
%% initial parameters 初始参数%
fa=9677; % 分析频率x
omega=2*pi*fa; % 角速度
c=340; k=omega/c; rho=1.225;pref=2e-5;
r=0.185;         % 半径
S=pi*r^2;        % 面积
zH = 0.45;       % Z coordinate of the array (of the measurements)
NumMic = 12;     % number of microphones
NumSM= 30;        % number of sequential measurements
Nsnap = 100;       % number of snapshots
load('Kappa.mat');
Kappa = Kappa/r; 
Kappa=Kappa(:,1); % 只考虑周向模态

%% The propagation mode numbers 根据截至频率，计算出可传播模态%
mode_prop2=propagated_models(k,Kappa);  % 可传播模态
[row,col] = size(mode_prop2); % 可传播模态数量

%% 扬声器阵列位置%
% load('speaker.mat');%传声器安装位置，要求柱坐标
% % figure
% % xs = speaker(:,1).*cos(speaker(:,2));ys=speaker(:,1).*sin(speaker(:,2));
% % scatter3(xs,ys,speaker(:,3),'LineWidth',1.5);

%% 传声器阵列%%%%% （柱坐标）
 mic_loc = zeros(NumSM*NumMic,3);
        XM=0.185*ones(NumMic*NumSM,1);
        YM=zeros(NumMic*NumSM,1);
        for  j=1:NumSM
        theta1=[0:30:330]';
         YM((1:NumMic)+(j-1)*NumMic,1)=theta1+(j-1)*1;
        end
         YM=YM/180*pi;
        ZM1=0.34*ones(NumMic,1);
        for  j=1:NumSM
            ZM((1:NumMic)+(j-1)*NumMic,1)=ZM1;
        end
        mic_loc=[XM,YM,ZM];
 
        figure
        col = ['ko';'ro';'go';'co';'mo';'yo'];
        col = [col;'kd';'rd';'gd';'cd';'md';'yd'];
        col = [col;'kx';'rx';'gx';'cx';'mx';'yx'];
        col = [col;'k^';'r^';'g^';'c^';'m^';'y^'];
        col = [col;'k>';'r>';'g>';'c>';'m>';'y>'];
        col = [col;'k<';'r<';'g<';'c<';'m<';'y<'];
        col = [col;'kv';'rv';'gv';'cv';'mv';'yv'];
        col = [col;'k*';'r*';'g*';'c*';'m*';'y*'];
         xs = mic_loc(:,1).*cos(mic_loc(:,2));ys=mic_loc(:,1).*sin(mic_loc(:,2));
        for i = 1:NumSM
      scatter3(xs((1:NumMic)+(i-1)*NumMic),ys((1:NumMic)+(i-1)*NumMic),mic_loc((1:NumMic)+(i-1)*NumMic,3),[col(1+mod(i-1,length(col)),:)],'LineWidth',1.5);
         hold on;
           end

        
%% 令目标模态的模态系数为一个定值%
      m_p = 2; n_p = 0;%%%%%%需要产生的模态阶数
      ind = ismember(mode_prop2,[m_p,n_p],'rows')';
      posi = find(ind==1);
      AA=zeros(1,row);AA(posi)=1;%%%%%目标模态系数为1，其他模态系数为0%%
      
%% 计算麦克风处的声压值%
       [P1]=presure1(mic_loc,mode_prop2,Kappa,k,r,AA);
        P1 =P1';
        a = normrnd(1,1,1,Nsnap);
        Pdirect=P1*a;
        SNR=30;
        Pnoise = AddNoise(Pdirect,SNR,NumSM,NumMic,Nsnap);     %         对声压信号添加噪声
        P=Pnoise;                              
        Spp1 = P*P'/Nsnap;
        SppOrg = (Spp1 + Spp1')./2;             % 相当于用 NumSM*Single 个麦克风进行测量而获得的矩阵(加噪声) 

          figure;
          imagesc(abs(SppOrg));colormap(flipud(hot));
          colorbar;
          axis square;xlabel('Microphone Index');ylabel('Microphone Index');
          set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);       

        %% Simulation Signal Definition of Sequential/Non-synchronous Measurement)
        SppSmOrg = zeros(NumSM*NumMic);
        for i = 1:NumSM
             SppSmOrg((1:NumMic)+(i-1)*NumMic,(1:NumMic)+(i-1)*NumMic) = P((1:NumMic)+(i-1)*NumMic,:)*P((1:NumMic)+(i-1)*NumMic,:)'/Nsnap;
        end
        
%% 构建观测矩阵%%%%%%%%
         [G]=matrix_G_trial(mode_prop2,Kappa,k,r,mic_loc);
         [ar,ac]=find(G==Inf);
         G(ar,ac)=0.5*(G(ar,ac-1)+G(ar,ac+1));
         cond(G)
 %%  获取非同步测量空间基函数      
        k = 2*pi*fa/c;
          K_p =fix( NumSM.^0.5.*NumMic);
%              K_p =68;
        [U,S,V] = svd(G);
        Phi_basis  = U(:,1:K_p);
        psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
        D_measured = SppSmOrg;          % measured matrix
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

        D_measured = SppSmOrg;          % measured matrix
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

%% 互谱矩阵重构误差求解
        FitErr = 100*norm(SppSmFit - SppOrg,'fro')/norm(SppOrg,'fro');
%         FitErr = 100*norm(SppSmFit - SppOrg_3,'fro')/norm(SppOrg_3,'fro');
        disp(['自-互谱矩阵重建误差1：',num2str(FitErr),'%']);
        
 %% 模态识别声压列向量
% Spp=SppSmOrg; 
 Spp=R_matrix_1;
P_amplitude = sqrt(diag(Spp)); % 获取幅值
P_phase = angle(Spp(:,1)/Spp(1,1));  % 获取相位角
% P_phase = angle(Spp(:,1)/Spporg_2(1,1)); 
P_complex = P_amplitude.*exp(1i*P_phase);
        FitErr = 100*norm(P_complex - P,'fro')/norm(P,'fro');
        disp(['声压列向量重建误差1：',num2str(FitErr),'%']);
P=P_complex;
 figure;
imagesc(abs(Spp));colormap(flipud(hot));%caxis([0,2500]);
colorbar;
axis square;xlabel('Microphone Index');ylabel('Microphone Index');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);  


%% 方法:1：最小二乘法%%%%
        q_re1=(G'*G)^-1*G'*P;
        
figure
abs_q1=20*log10(abs(q_re1)/pref);%/(2*10-5)
abs_q1(find(abs_q1<0))=0;
bar(abs_q1)
b =bar(abs_q1,'FaceColor',[.8 .85 .95],'LineWidth',1);
b.FaceColor = 'flat';
b.CData(21,:) = [0 0 0.54];
xlabel('Modal order');ylabel('Modal coefficient（dB）');
xlim([0 size(mode_prop2,1)+1])
colormap(hot)
ylim([0 max(abs_q1)+15])
set(gca,'YTick',[0:35:140]) 
% title('最小二乘法')
set(gca,'XTick',1:row);
   set(gcf,'position',[50 400 1400 300]);
 
   %% 频率为12084Hz时周向模态的排列
% set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)','(15,0)','(16,0)','(17,0)','(18,0)','(19,0)','(20,0)','(21,0)','(22,0)','(23,0)','(24,0)','(25,0)','(26,0)','(27,0)','(28,0)','(29,0)','(30,0)','(31,0)','(32,0)',...
%  '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
% '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)','(-15,0)','(-16,0)','(-17,0)','(-18,0)','(-19,0)',...
% '(-20,0)','(-21,0)','(-22,0)','(-23,0)','(-24,0)','(-25,0)','(-26,0)','(-27,0)','(-28,0)','(-29,0)','(-30,0)','(-31,0)','(-32,0)'})

%% 频率为9667Hz时周向模态排列
  set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)','(15,0)','(16,0)','(17,0)','(18,0)','(19,0)','(20,0)','(21,0)','(22,0)','(23,0)','(24,0)','(25,0)','(26,0)','(27,0)','(28,0)','(29,0)',...
   '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
 '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)','(-15,0)','(-16,0)','(-17,0)','(-18,0)','(-19,0)',...
'(-20,0)','(-21,0)','(-22,0)','(-23,0)','(-24,0)','(-25,0)','(-26,0)','(-27,0)','(-28,0)','(-29,0)'});
        %% 频率为4833Hz时周向模态排列
%           set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)',...
%            '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
%          '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)'});
%%  频率为2900Hz时周向模态的排列
% set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(-1,0)','(-2,0)','(-3,0)',...
%    '(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)'});
set(gca,'FontSize',14)
%% 方法:2：GMC%%%%
% lamda=0.01;
% zeta=10^(-0.3);
% rho=200*lamda*zeta;
% gamma=1.1;
% sc=1*10^-4;
% mlter=1000;
% [para_val,q_re2,time]=ADMM_SRLASSO(G,P,lamda,rho,gamma,zeta, sc,mlter);
% 
% figure
% abs_q2=20*log10(abs(q_re2)/pref);%/(2*10-5)
% abs_q2(find(abs_q2<0))=0;
% bar(abs_q2)
% xlabel('Modal order');ylabel('Modal coefficient （dB）');
% xlim([0 size(mode_prop2,1)+1])
% colormap(hot)
%  ylim([max(abs_q2)-15 max(abs_q2)+15])
% title('GMC')
% set(gca,'XTick',1:row);
% % set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)'...
% %     '(0,1)','(1,1)','(2,1)','(3,1)','(4,1)','(5,1)','(6,1)','(0,2)','(1,2)','(2,2)','(3,2)','(0,3)','(1,3)','(-1,0)','(-2,0)','(-3,0)',...
% %     '(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)','(-9,0)','(-10,0)','(-1,1)','(-2,1)','(-3,1)','(-4,1)','(-5,1)','(-6,1)','(-1,2)','(-2,2)','(-3,2)','(-1,3)'})
% set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(-1,0)','(-2,0)','(-3,0)',...
%    '(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)','(-9,0)','(-10,0)'});

%% 方法3：贝叶斯压缩感知
         PHI=G;
         eta=1e-8;
         adaptive=0;
         optimal=1;
         scale=0.1;
         t=P;
         sigma2 =std(t)^2/1e2;
         q_re6= zeros(row,1);
         [weights,used,sigma2,errbars] = BCS_fast_rvm(PHI,t,sigma2,eta,adaptive,optimal,scale)
         q_re6(used)= weights;
         
figure
abs_q6=20*log10(abs(q_re6)/pref);%/(2*10-5)
abs_q6(find(abs_q6<0))=0;
bar(abs_q6)
b =bar(abs_q6,'FaceColor',[.8 .85 .95],'LineWidth',1);
b.FaceColor = 'flat';
b.CData(21,:) = [0 0 0.54];
xlabel('Modal order');ylabel('Modal coefficient（dB）');
xlim([0 size(mode_prop2,1)+1])
colormap(hot)
ylim([0 max(abs_q1)+15])
set(gca,'YTick',[0:35:140]) 
% title('贝叶斯压缩感知') 
 set(gcf,'position',[50 400 1400 300]);
set(gca,'XTick',1:row);

   %% 频率为12084Hz时周向模态的排列
% set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)','(15,0)','(16,0)','(17,0)','(18,0)','(19,0)','(20,0)','(21,0)','(22,0)','(23,0)','(24,0)','(25,0)','(26,0)','(27,0)','(28,0)','(29,0)','(30,0)','(31,0)','(32,0)',...
%  '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
% '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)','(-15,0)','(-16,0)','(-17,0)','(-18,0)','(-19,0)',...
% '(-20,0)','(-21,0)','(-22,0)','(-23,0)','(-24,0)','(-25,0)','(-26,0)','(-27,0)','(-28,0)','(-29,0)','(-30,0)','(-31,0)','(-32,0)'})

%% 频率为9667Hz时周向模态排列
  set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)','(15,0)','(16,0)','(17,0)','(18,0)','(19,0)','(20,0)','(21,0)','(22,0)','(23,0)','(24,0)','(25,0)','(26,0)','(27,0)','(28,0)','(29,0)',...
   '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
 '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)','(-15,0)','(-16,0)','(-17,0)','(-18,0)','(-19,0)',...
'(-20,0)','(-21,0)','(-22,0)','(-23,0)','(-24,0)','(-25,0)','(-26,0)','(-27,0)','(-28,0)','(-29,0)'});
        %% 频率为4833Hz时周向模态排列
%           set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)',...
%            '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
%          '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)'});
%% 频率为2900Hz时周向模态的排列
% set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(-1,0)','(-2,0)','(-3,0)',...
%    '(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)'});
set(gca,'FontSize',14)






























































