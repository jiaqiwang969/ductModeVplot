% Aim: 利用试验数据生成互谱矩阵
% % 输出整个频率

clc;
clear;
close all;

%% add subfunction
addpath(genpath('.'));
addpath(genpath('../'));
chemin = '../database/03-16mic';

%% add Basic parameters
algorithm=3;

zH = 0;             % 测试距离
zRef=zH-0.08;       % 参考传声器到测试面的轴向距离
NumMic = 12;        % 传声器的数量
NumSM  = 30;        % 测量的次数
Radius = 0.2;       % 管道半径
Area = pi*Radius^2; % 管口面积
Fs = 25600;         % 采样频率
time = 5;           % 采样时间

rotor_speed=8000;   % 轴转速信息

%w=2*pi*rotor_speed/60*29/343*Radius;

%% data processing
L_signal = Fs*time;             %信号长度
L_seg = round(L_signal/100);    %确定对信号处理划窗长度
Wind =  hamming(L_seg);         %确定对数据进行汉宁窗处理
Noverlap = round(L_seg/2);      %确定信号划窗重叠率
Nfft = 2^(ceil(log2(L_seg))+1); %确定分析频率
round=7;                        %分段,生成多个block，每个block 为round
name='m8_3400Hz_Rotate';

%% 方法01:CC算法
Ind = [1:NumSM];   %设定循环次数
data_all=[];
Num_file = Ind ;
for i_file =Num_file
    eval(['load ''',chemin,'/',name,'-No-',num2str(i_file),'.mat''']);       %读取数据
    data_all=[data_all Data(:,1:13)];
    Tdata{i_file}=Data(:,1:13);
    [temp_ref,freq] = cpsd(Tdata{i_file}(:,13),Tdata{i_file}(:,13),Wind,Noverlap,Nfft,Fs);
    temp_ref = sqrt(temp_ref);
    for k=1:NumMic
        [temp,freq] = cpsd(Tdata{i_file}(:,k),Tdata{i_file}(:,13),Wind,Noverlap,Nfft,Fs);
        CC(:, i_file+NumSM*(k-1)) = temp./temp_ref;                                                 %360个测点的数据矩阵
    end
end

NumMic_enlarge=NumSM*NumMic;                                                                             %传声器可测量的模态总数
mode=[-60:60];
for k=1:length(mode)
    a_mf(:,k)=1/NumMic_enlarge*CC(:,1:NumMic_enlarge)*exp(mode(k)*1i*2*pi*(0:NumMic_enlarge-1)/NumMic_enlarge).';    %求解模态系数
end
GAMMA =abs(a_mf);


%% 方法02:合成孔径算法
%% 是否将参考传声器植入到互功率谱矩阵中
for if_consider_ref=[0 1] %with method1: reference sensor is same to others
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

    df =freq(2) - freq(1);
    f0=[30:10:3400*2.2];
    Freq_slice=f0/(rotor_speed/60*29);
    for k=1:length(Freq_slice)
        xunhao_around=floor(f0(k)/df)+[-3:3];
        xuhao=xunhao_around(1)+find(amf(xunhao_around)==max(amf(xunhao_around)))-1;
        CSM{k}=[];
        for i_file=1:NumSM
            CSM{k} = blkdiag(CSM{k},reshape(CC3(xuhao,:,:,i_file),12+if_consider_ref,12+if_consider_ref));
        end
        % 对CC进行处理，形成真正对ref参考矩阵
        if if_consider_ref==1
            % Step01: 裁掉传声器所在对行
            CC_1=CSM{k};
            CC_1(13*[1:30],:)=[];
            CC_1(:,13*[1:30])=[];
            % Step02: 补上参考传声器的影响
            CC_2=CSM{k};CC_2(13*[1:30],:)=[]; temp1=sum(CC_2(:,13*[1:30]),2);
            CC_3=CSM{k};CC_3(:,13*[1:30])=[]; temp2=sum(CC_3(13*[1:30],:),1);
            CC_diag=diag(CSM{k});temp3=mean(CC_diag(13*[1:30]));
            % Step03: 重新组装
            CSM{k}=[[temp3 temp2];[temp1 CC_1]];
        end
    end

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

    % 根据截至频率，计算出可传播模态
    parfor k=1:length(Freq_slice)

        %% 该处替换为wjq写的green函数
        % 输入：w,Radius,mic_loc
        % 输出：G:361*[-mn,+mn] 频率域
        k
        m = [-60:60];   % 周向模态限制范围，自动删选
        n = [1];        % 径向模态限制范围
        M =  0;       % 管道流速
        [G,index_mn{k}]=matrix_G_basis(f0(k),Radius,M,mic_loc,m,n);
        %cond(G)


        %% 非同步测量空间基函数的确定
        [U,S,V] = svd(G);
        Phi_basis  = U;
        psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
        D_measured = CSM{k};                % measured matrix

        %% 非同步测量算法代码
        %algorithm=1,2,3对应选择不同算法
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
                SC = 0.05;                      % stopping criteria
                mIter = 34;                      % maximum iteation
                gama =5;       %relaxation parameter
                alpha = 28.5e-3; %regularization parameter
                mu=29/n;       %penalty parameter
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
                tic
                [R_matrix_1,para_AP]  = CP(Hard_Th,D_measured,psi_B,max_it,SC,meanORG,varORG);
                toc;

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
        q_re1{k}=(G'*G)^-1*G'*P;   %03-30
    end

    if if_consider_ref==0
        q_re1_noref=q_re1;
    else
        q_re1_ref=q_re1;
    end

end

Mode1=[[-60:60].'];
Qr_noref=zeros(length(Mode1),length(index_mn));
Qr_ref=zeros(length(Mode1),length(index_mn));

for k=1:length(index_mn)
    tp1=size(index_mn{1,k},1);
    Qr_noref(61-(tp1-1)/2:61+(tp1+1)/2-1,k)=q_re1_noref{k};
    Qr_ref(61-(tp1-1)/2:61+(tp1+1)/2-1,k)=q_re1_ref{k};

end

Gamma_noref=flip(abs(Qr_noref).',2);
Gamma_ref=flip(abs(Qr_ref).',2);


figure
subplot(2,1,1)
bar(mode,GAMMA(find(freq==3400),:))
xlim([-10,10])

% subplot(2,1,2)
% bar([-60:60],Gamma_ref(find(f0==3400),:))
% xlim([-10,10])




h=figure;
subplot(3,2,1)
bar(mode,GAMMA(find(freq==3400),:))
%ylim([0,3400*2.2])
axis xy
grid minor
% colorbar
xlim([-11,11])
title('Tonal')
subplot(3,2,2)
imagesc(mode,freq,GAMMA./sum(GAMMA,2))
ylim([0,3400*2.2])
axis xy
colorbar
xlim([-25,25])
title('Normalized')
subplot(3,2,3)
bar(Mode1,Gamma_noref(find(f0==3400),:))
% imagesc(Mode1,f0,Gamma_noref)
% ylim([0,3400*2.2])
axis xy
grid minor
% colorbar
xlim([-11,11])
subplot(3,2,4)
imagesc(Mode1,f0,Gamma_noref./sum(Gamma_noref,2))
axis xy
colorbar
xlim([-25,25])
ylim([0,3400*2.2])
subplot(3,2,5)
bar(Mode1,Gamma_ref(find(f0==3400),:))
% imagesc(Mode1,f0,Gamma_ref)
%ylim([0,3400*2.2])
axis xy
grid minor
% colorbar
xlim([-11,11])
subplot(3,2,6)
imagesc(Mode1,f0,Gamma_ref./sum(Gamma_ref,2))
axis xy
colorbar
ylim([0,3400*2.2])
xlim([-25,25])


saveas(h,[date,name,'-',num2str(algorithm),'.fig'])
saveas(h,[date,name,'-',num2str(algorithm),'.png'])




