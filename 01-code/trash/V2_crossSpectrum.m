clc;
clear;
close all;

%% add subfunction
addpath(genpath('.'));
chemin = '../database/01-rotateMic';

%% add Basic parameters
a=0.185;           % 管道半径
zH = 0.4;          % 阵列的Z坐标
NumMic = 12;           % 传声器的数量
NumSM= 30;         % 非同步测量的次数
Fs = 102400 ;      % 采样频率
time=5;            % 采样时长
Nw =  1024*5;      % length of snapshot，resolution =   Fs/Nw；
nfft =2*(Nw);      % nfft = 2^ (Nw);
F = Fs*(0:nfft-1)/nfft;   % 频率画点
w = hanning(Nw);    % window function
w = w/norm(w);      % w = w*sqrt(2/nfft); % calibration
rotor_speed = 12000;        %轴转速信息
F0=rotor_speed/60*29;    %Frequency
[temp,indice_f] = min(abs(F-rotor_speed));       %找到最接近波束形成频率f的频率所在列数
Nsnap = 100;        % Nsnap = floor(2*((Fs*10))/Nw-1);%快照数量

%% 构建互谱矩阵
pM_temp = zeros(NumMic,Nsnap);
pM = [];
p=zeros(Fs*time,NumMic);
Ind = [1:NumSM];
Num_file = Ind ;

for k =Num_file
    eval(['load ''',chemin,'/','RotaryTest-12000-Rotate-No-',num2str(k),'.mat''']);
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
axis equal



%% 构建模态系数矩阵%%%%%%%%
nk_enlarge=NumSM*NumMic;
m=-nk_enlarge/2:nk_enlarge/2;
for k=1:length(m)
    G(:,k)=exp(m(k)*-1i*2*pi*(1:nk_enlarge)/nk_enlarge).';
end

%% 非同步测量空间基函数的确定
K_p =fix( NumSM.^0.5.*NumMic);
%          K_p =120;
%         [U,S,V] = svd(G);
Phi_basis  = G(:,1:K_p);
psi_B = Phi_basis*pinv(Phi_basis'*Phi_basis)*Phi_basis';
D_measured = MM;          % measured matrix


%% 非同步测量算法
% ADMM算法
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
bar([-nk_enlarge/2:nk_enlarge/2],abs_q1)
colormap(hot)
ylim([0 max(abs_q1)+15])
set(gca,'YTick',[0:35:140])
% title('最小二乘法')
set(gcf,'position',[50 400 800 300]);






