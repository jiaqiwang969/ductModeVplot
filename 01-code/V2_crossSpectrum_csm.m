% 利用fft变换得到的CSM，
% 参考https://github.com/Jiaqi-knight/RotaryCasingWallMeasure-Testing/blob/3042d9ea599d2df46df8bcbb9489c548ba0d5d71/CSM.m
% 输出单个频率

% clc;
% clear;
% close all;

chemin = '../database/01-rotateMic';
%% 基础参数设置
r=0.185;           % 管道半径
rhoAir = 1.2;      % 空气密度
c = 340;           % 声速
zH = 0.34;         % 阵列的Z坐标
NumMic = 12;       % 传声器的数量
NumSM= 30;         % 非同步测量的次数
f0 = 12000/60*29*1;         % 分析频率
Fs = 102400 ;      % 采样频率
time=5;            % 采样时长
Nw =  1024*5;      % length of snapshot，resolution =   Fs/Nw；
nfft =2*(Nw);      % nfft = 2^ (Nw);
F = Fs*(0:nfft-1)/nfft;   % 频率画点
w = hanning(Nw);   % window function
w = w/norm(w);     % w = w*sqrt(2/nfft); % calibration
[temp,indice_f] = min(abs(F-f0));       %找到最接近波束形成频率f的频率所在列数
Nsnap = 100;       % Nsnap = floor(2*((Fs*10))/Nw-1);%快照数量


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

