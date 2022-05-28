clc;
clear;
close all;
addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\'));
addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\simulation\'));
addpath(genpath('E:\模态识别\模态识别代码\mode_detect\mode_detect\simulation\subprogram\'));
chemin = 'E:\模态识别\基于非同步测量的压气机管道模态识别数据\实验20-2019-11-16旋转机闸 测试';
%% 基础参数设置
        r=0.185;           % 管道半径
        rhoAir = 1.2;      % 空气密度
        c = 340;           % 声速
        zH = 0.34;         % 阵列的Z坐标
        NumMic = 12;       % 传声器的数量
        NumSM= 2;         % 非同步测量的次数
        f0 = 2900*2;         % 分析频率
        Fs = 204800 ;      % 采样频率
        time=2;            % 采样时长
        Nw =  1024*5;        % length of snapshot，resolution =   Fs/Nw；
        nfft =2*(Nw);          % nfft = 2^ (Nw);
        F = Fs*(0:nfft-1)/nfft;   % 频率画点
        w = hanning(Nw);    % window function
         w = w/norm(w);      % w = w*sqrt(2/nfft); % calibration
        [temp,indice_f] = min(abs(F-f0));       %找到最接近波束形成频率f的频率(0,所在列数) -1000：10：48990
        Nsnap = 100;        % Nsnap = floor(2*((Fs*10))/Nw-1);%快照数量
        
%% 传声器阵列%%%%% （柱坐标）
        mic_loc = zeros(NumSM*NumMic,3);
        XM=r*ones(NumMic*NumSM,1);
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

        %% 构建互谱矩阵        
        pM_temp = zeros(NumMic,Nsnap);
        pM = [];
        p=zeros(Fs*time,NumMic);
        Ind = [1:NumSM];   
         Num_file = Ind ; 
        
       for k =Num_file
    eval(['load ''',chemin,'\','RotaryTest-2-6000-Rotate-No-',num2str(k),'.mat''']);
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
%  load('E:\非同步测量\ALM_ADMM非同步测量最终代码\22.04.24_Data\MM_2900.mat');
%% 选择管道内可传播的模态%%%%%%%%
omega=2*pi*f0; % 角速度
k=omega/c; % 波数
load('Kappa.mat');
Kappa = Kappa/r; 
Kappa=Kappa(:,1); % 只考虑周向模态

%% 根据截至频率，计算出可传播模态%
          mode_prop2=propagated_models(k,Kappa);  % 可传播模态
          [row,col] = size(mode_prop2);           % 可传播模态数量
         [G]=matrix_G_trial(mode_prop2,Kappa,k,r,mic_loc);
%          [ar,ac]=find(G==Inf);
%          G(ar,ac)=0.5*(G(ar,ac-1)+G(ar,ac+1));
         cond(G)
         
 %% 非同步测量空间基函数的确定         
        K_p =fix( NumSM.^0.5.*NumMic);
%          K_p =120;
        [U,S,V] = svd(G);
        Phi_basis  = G(:,1:35);
        psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
        D_measured = MM;          % measured matrix
        
        
%% 非同步测量算法
        
% ADMM
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
        
        %% 由互谱矩阵获取声压列向量
%         Spp=MM;
        Spp=R_matrix_1;
        P_amplitude = sqrt(diag(Spp));
        P_phase = angle(Spp(:,1)/Spp(1,1));
        P_complex = P_amplitude.*exp(1i*P_phase);
        P=P_complex;

        %% 模态识别方法:1：最小二乘法%%%%
        q_re1=(G'*G)^-1*G'*P;
         pref=2e-5;
         
figure
        abs_q1=20*log10(abs(q_re1)/pref);%/(2*10-5)
        abs_q1(find(abs_q1<0))=0;
        bar(abs_q1)
        b =bar(abs_q1,'FaceColor',[0 0 0.54],'LineWidth',1);
        b.FaceColor = 'flat';
        xlabel('Modal order');ylabel('Modal coefficient（dB）');
        xlim([0 size(mode_prop2,1)+1])
        colormap(hot)
        ylim([0 max(abs_q1)+15])
        set(gca,'YTick',[0:35:140]) 
        % title('最小二乘法')
         set(gca,'XTick',1:row);
        set(gcf,'position',[50 400 800 300]);
        %%  频率为2900Hz时周向模态的排列
%         set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(-1,0)','(-2,0)','(-3,0)',...
%            '(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)'});
        %% 频率为3866Hz时周向模态排列
%           set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)',...
%            '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
%          '(-9,0)','(-10,0)','(-11,0)'});
        %% 频率为4833Hz时周向模态排列
%           set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)',...
%            '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
%          '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)'});
        %% 频率为9667Hz时周向模态排列
   set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)','(15,0)','(16,0)','(17,0)',...
           '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
         '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)','(-15,0)','(-16,0)','(-17,0)'});
        set(gca,'FontSize',14)      

%% 方法2：贝叶斯压缩感知
%          PHI=G;
%          eta=1e-8;
%          adaptive=0;
%          optimal=1;
%          scale=0.1;
%          t=P;
%          sigma2 = std(t)^2/1e2;
%          q_re6= zeros(row,1);
%          [weights,used,sigma2,errbars] = BCS_fast_rvm(PHI,t,sigma2,eta,adaptive,optimal,scale)
%          q_re6(used)= weights;
% 
% 
% figure
%         abs_q6=20*log10(abs(q_re6)/pref);%/(2*10-5)
%         abs_q6(find(abs_q6<0))=0;
%         bar(abs_q6)
%         b =bar(abs_q6,'FaceColor',[0 0 0.54],'LineWidth',1);
%         b.FaceColor = 'flat';
%         xlabel('Modal order');ylabel('Modal coefficient（dB）');
%         xlim([0 size(mode_prop2,1)+1])
%         colormap(hot)
%         ylim([0 max(abs_q1)+15])
%         set(gca,'YTick',[0:35:140]) 
%         % title('贝叶斯压缩感知') 
%          set(gcf,'position',[50 400 800 300]);
%         set(gca,'XTick',1:row);
% 
% 
%         %% 频率为9667Hz时周向模态排列
%         %   set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)','(15,0)','(16,0)','(17,0)','(18,0)','(19,0)','(20,0)','(21,0)','(22,0)','(23,0)','(24,0)','(25,0)','(26,0)','(27,0)','(28,0)','(29,0)',...
%         %    '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
%         %  '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)','(-15,0)','(-16,0)','(-17,0)','(-18,0)','(-19,0)',...
%         % '(-20,0)','(-21,0)','(-22,0)','(-23,0)','(-24,0)','(-25,0)','(-26,0)','(-27,0)','(-28,0)','(-29,0)'});
%         %% 频率为4833Hz时周向模态排列
%           set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)','(12,0)','(13,0)','(14,0)',...
%            '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
%          '(-9,0)','(-10,0)','(-11,0)','(-12,0)','(-13,0)','(-14,0)'});
%         %% 频率为3866Hz时周向模态排列
%         %   set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(9,0)','(10,0)','(11,0)',...
%         %    '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)',...
%         %  '(-9,0)','(-10,0)','(-11,0)'});
%         %% 频率为2900Hz时周向模态的排列
% %         set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)','(-1,0)','(-2,0)','(-3,0)',...
% %            '(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)'});
%         set(gca,'FontSize',14)

 