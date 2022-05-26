clc
clear

Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 150000;           % Length of signal
t = (0:L-1)*T;        % Time vector
N=10;
time=L/Fs;
for k=1:10
    X(:,k) = sin(2*pi*50*t+2*pi/N*k) ;
end
figure
plot(1000*t(1:50),X(1:50,:))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')
Y = fft(X);
P2 = real(Y/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
f = Fs*(0:(L/2))/L;
figure
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% 验证 CPSD
L_signal = Fs*time;             %信号长度
L_seg = round(L_signal/10);     %10截断
Wind = hamming(L_seg);          %确定对数据进行汉宁窗处理
Noverlap = round(L_seg/2);      %确定信号划窗重叠率
Nfft = 2^(ceil(log2(L_seg))+1); %确定分析频率
CC1=[];
T1=  kron(ones(1,N), X  );
T2=  kron(X, ones(1,N)  );
[temp,freq]=cpsd(T1,T2,Wind,Noverlap,Nfft,Fs);
CC1=[CC1 temp];
CC3=reshape(CC1,length(freq),10,10);
f0=50; %观察频率/Hz
df=freq(2)-freq(1);
xunhao_around=floor(f0/df)+[-3:3];
amf=CC1(:,1);
xuhao=xunhao_around(1)+find(amf(xunhao_around)==max(amf(xunhao_around)))-1;
CC=permute(CC3(xuhao,:,:),[2,3,1]);
figure;
subplot(1,2,1)
imagesc(abs(CC))  %? 幅值如何修正？ /5.2174
colorbar;axis equal
title('cpsd')

%% 验证 CSM
Nw =  1000*1;      % length of snapshot，resolution =   Fs/Nw；
nfft =2*(Nw);      % nfft = 2^ (Nw);
F = Fs*(0:nfft-1)/nfft;   % 频率画点
w = hanning(Nw);   % window function
w = w/norm(w);     % w = w*sqrt(2/nfft); % calibration
f0=50;
[temp,indice_f] = min(abs(F-f0));       %找到最接近波束形成频率f的频率所在列数
Nsnap = 10;        % Nsnap = floor(2*((Fs*10))/Nw-1);%快照数量

pM = [];
p=zeros(Fs*time,10);
for index_M = 1:10
    p(:,index_M)=X(:,index_M);
    for m = 1:Nsnap
        Pf = conj(fft(p((m-1)*Nw/2+(1:Nw),index_M).*w,nfft));
        pM_temp(index_M,m) = Pf(indice_f);
    end
end
pM = pM_temp;
Cov_pM = pM*pM'/Nsnap;
Cov_pM = (squeeze(Cov_pM) + squeeze(Cov_pM)')/2;%squeeze函数用于删除矩阵中的单一维
MM=Cov_pM;
subplot(1,2,2)
imagesc(abs(MM));
colorbar;axis equal
title('csm')





