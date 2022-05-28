%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  ѭ��ά���˲���
%               �����������¶ԣ�����ѭ��ƽ�ȣ���������������ȡ
% ��Ҫ�ο����ף�
% (1)Liang Yu, etc, Extraction and imaging of aerodynamically generated sound field of rotor blades 
% in the wind tunnel test, MSSP2018
% (2)Liang Yu, etc, Reconstruction of cyclostationary sound source based on a back-propagating
% cyclic wiener filter, JSV2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
clear
close all
% �ļ�·��
SampeFreq = 51200 ;  % ����Ƶ��
AnaChannelNums = 1:135;  % ������10ͨ��
SampleTime = 1; % ����ʱ��
% ����ԭʼ��������
[Datafilename,mypath1]=uigetfile('D:\Aѧϰ����\����������\��������3\Bo105\��ͣ\*.*','���뽰Ҷ��˷��������'); 
S=strcat(mypath1,Datafilename);
fid = fopen(S,'r','b');
Binaryfs=fread(fid,1,'single');
Binarynummic=fread(fid,1,'single');
CLenth=fread(fid,1,'single');
NLenth=fread(fid,1,'single');
RawData=zeros(CLenth*NLenth,Binarynummic);
tic
for NNi=1:NLenth
    RawData(((NNi-1)*CLenth+1):(NNi*CLenth),:)=fread(fid,[CLenth Binarynummic ],'single');
end
toc
fclose(fid);
disp(['����������ݶ���!'])
%% ��������
y = RawData(1:SampleTime*SampeFreq,AnaChannelNums); % ��ԭʼ�����ж�ȡ��������
L = size(y,1);  % �������ݵĳ���
Fs = SampeFreq;
% ��ʾ�������ݺͻ�����Ƶ�׷���
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
[S,f] = pwelch(y(:,end),Nw,Nv,Nfft,Fs);
figure
subplot(211)
plot([0:1/Fs:1-1/Fs],y(1:Fs,1));
xlabel('ʱ�� (��)'),ylabel('����'),
title('��ʾ��һ���ڵĲ�������');
subplot(212)
plot(f(1:6000),10*log10(S(1:6000))) % ֻ��ʾ��6000Hz
xlabel('Ƶ�� (����)'),ylabel('dB'),
title(['�ź�Ƶ�� ; Ƶ�ʷֱ��� = ',num2str(diff(f(1:2))),'Hz'])

% ѭ��Ƶ������Է���
alpha_max = 300; % ����ѭ��Ƶ��
Nw = 2^8; % ���ڳ���
opt.coh = 1;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);
figure
plot(alpha(2:end),sum(abs(S(:,2:end)))),
xlabel('ѭ��Ƶ�� (Hz)')
title('ѭ��Ƶ����ʾ');

%%  ʱ��ͬ��ƽ����ȥ��һ��г����
% [Datafilename,mypath1]=uigetfile('D:\Aѧϰ����\����������\��������3\Bo105\ƽ��\*.*','���뽰Ҷ��λ����'); 
% S=strcat(mypath1,Datafilename);
% fid = fopen(S,'r','b');
% Nona=fread(fid,1,'single');
% Timenummic=fread(fid,1,'single');
% CLenthTime=fread(fid,1,'single');
% NLenthTime=fread(fid,1,'single');
% BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
% fclose(fid);
% disp(['��ɷ�λ�����ݶ���!'])
% time_blade = BladePositionSig;
% 
%  ͬ��ƽ����ȡһ��ѭ��ƽ���ź�
% [~, syn_c_num] = size(time_blade);
% ave_syn_y = zeros(size(y));
% y_residual = zeros(size(y));
ave_syn_y = zeros(49698,135);
y_residual = zeros(49698,135);

for i = 1:135
    time_blade = 136;
    Period_time = 1/136;
    Period_sample = round(Fs*Period_time*6);
    Period_number = fix(L/Period_sample);
    y_Period = y(1:Period_number*Period_sample,i);
    y_Period_matrix = reshape(y_Period,Period_sample,Period_number);
    My_Period = mean(y_Period_matrix,2);
    ave_syn_y(:,i) = repmat(My_Period,Period_number,1);
    y_residual(:,i) = y_Period - ave_syn_y(:,i);
end
%%  �������˲�����ȥ��һ��г���źţ�
% % ��ȡת���ļ�
% [Datafilename,mypath1]=uigetfile('D:\Aѧϰ����\����������\��������3\Bo105\ƽ��\*.*','���뽰Ҷ��λ����'); 
% S=strcat(mypath1,Datafilename);
% fid = fopen(S,'r','b');
% Nona=fread(fid,1,'single');
% Timenummic=fread(fid,1,'single');
% CLenthTime=fread(fid,1,'single');
% NLenthTime=fread(fid,1,'single');
% BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
% fclose(fid);
% disp(['��ɷ�λ�����ݶ���!'])
% time_blade = BladePositionSig;
%     
% % ����ת����Ϣͬ��һ��˲ʱƵ��
% f_rot=[]; % ˲ʱƵ��
% for i=1:length(time_blade)
%     f_i=1/time_blade(i); % ��Ƶ��1/time_blade(i) ��ҶƵ��4�ǽ�Ҷ�ĸ���
%     f_rot=[f_rot f_i*ones(1,round(time_blade(i)*Fs))];
% end
% f_rot = f_rot.';

% % �ź�׼��
% f_rot =34;
% fp = f_rot*ones(size(y,1),1); 
% % t = (1:size(y_xhat,2))/Fs; %��ɢʱ���ź�
% 
% % �������˲�������׼��
% ave_syn_y = zeros(size(y));
% y_residual = zeros(size(y));
% for i=1:135
%     r = 2000;   % ����Ȩ�� 2000
%     ford = 1;    % �����˲�������Ϊ1��Ҳ����ȡ2��
%     tol = 0.1;    % ����
%     maxit = 500;  % ����������
%     V_fp = [4*fp 2*4*fp 3*4*fp 4*4*fp 5*4*fp 6*4*fp 7*4*fp 8*4*fp 9*4*fp 10*4*fp 11*4*fp 12*4*fp 13*4*fp 14*4*fp 15*4*fp 16*4*fp 17*4*fp 18*4*fp 19*4*fp 20*4*fp 21*4*fp 22*4*fp 23*4*fp 24*4*fp 25*4*fp 27*4*fp 28*4*fp 29*4*fp 30*4*fp 31*4*fp 32*4*fp 33*4*fp 34*4*fp 35*4*fp 37*4*fp 38*4*fp 39*4*fp 40*4*fp  0*fp]; % 4 ��ҶƬ�ĸ���, 5 ��ϣ���������г������
%     [nt,nch] = size(V_fp);
%     [xm,bwm,Tm,fl,rr,it,rv,xr] = vkm(y(:,i)',V_fp,Fs,r*ones(nch,1),ford,tol,maxit);  
%     ave_syn_y(:,i)=sum(xr(:,1:39),2);
%     y_residual(:,i) = y(:,i)- ave_syn_y(:,i);    
% end

%% �źŷ��� 
% % FFTƵ�׷���
% frequency_disp = Fs/2;  %��ʾƵ�ʵ�**Hz
% f = Fs*(0:(L/2))/L;
% Y = fft(y_residual(:,1)); %FFT
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% figure
% plot(f,P1,'k','LineWidth',1.5) 
% title('FFT Amplitude Spectrum','fontweight','b')
% xlabel('f [Hz]')
% ylabel('Amplitude')    % ע�⣺����ϵ�λ[m/s^2]��[Pa]
% xlim([0 frequency_disp])
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% set(gcf,'position',[400 250 500 250])

% �����׷�����Welch��s power spectral density estimate��
Nw =Fs/20;   % window
Nv =  ceil(3/4*Nw); % overlap
Nfft = Nw;  % nfft
frequency_disp = Fs/2;
figure
[S,f] = pwelch(y_residual(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp/20),10*log10(S(1:frequency_disp/20)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

%%  ѭ��ά���˲�����ȡ����ѭ��ƽ���źţ�
L1 = size(y_residual,1);
Na = 10; % ѭ��Ƶ�ʵĽ���
orders = (1:Na/2); % ѭ��Ƶ�ʵĸ���
orders = [orders;-orders]; 
orders = orders(:); % �ܵ�ѭ��Ƶ�ʸ���
alpha = 34.2/Fs; % ��һ��ѭ��Ƶ�� /136.8
phi = 2*pi*alpha*(0:L1-1)'; % ��λ������
fi = alpha*ones(L1,1); % Ƶ�ʵ�����

% ����ѭ�����׾���
Nw = 2^8;
[Syya,Syaya] = CyclicSpecMat(y_residual(:,:),phi,fi,orders,Nw);
f = Syaya.f*Fs;

% ����ѭ��ά���˲���
Ns = 1; % ѭ��ƽ����Դ�ĸ���
[G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,1e-6);

% ����Gabor�任�����˲�
[xhat,cx] = CyclicSpatialFiltering(y_residual(:,:),G(:,:,:),phi,fi,orders);
% [L_xhat,~] = size(xhat);
% nhat = y_residual(1:L_xhat) - xhat;
% y_residual=y_residual(1:51136,:);
% nhat = y_residual - xhat;
%% ����ͼ��
% ǰһ��Ĳ�����������ȡ�����ݽ��жԱ�
figure;
subplot(311);
plot([0:1/Fs:0.1-1/Fs],y(1:0.1*Fs,end));
xlabel('ʱ�� (��)'),ylabel('����'),
title('��������');
subplot(312);
plot([0:1/Fs:0.1-1/Fs],y_residual(1:0.1*Fs,end));
xlabel('ʱ�� (��)'),ylabel('����'),
title('ȥ��г��ʣ���ź�');
subplot(313);
plot([0:1/Fs:0.1-1/Fs],xhat(1:0.1*Fs,end));
xlabel('ʱ�� (��)'),ylabel('����'),
title('��ȡ��ѭ��ƽ���ź�');
% subplot(313);
% plot([0:1/Fs:0.1-1/Fs],nhat(1:0.1*Fs,end));
% xlabel('ʱ�� (��)'),ylabel('����'),
% title('�в��ź�ʱ������');

% �˲����ѭ��Ƶ������Է���
alpha_max = 300; % ����ѭ��Ƶ��
Nw = 2^8; % ���ڳ���
opt.coh = 1;

figure,
[S,alpha,f,STFT,t,Nv] = Fast_SC(xhat(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end)))), hold on;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y_residual(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end))),'r');
% [S,alpha,f,STFT,t,Nv] = Fast_SC(nhat(:,end),Nw,alpha_max,Fs,opt);
% plot(alpha(2:end),sum(abs(S(:,2:end))),'y');
xlabel('ѭ��Ƶ�� (Hz)')
legend('��ȡ���ź�','�������ź�')
% legend('��ȡ���ź�','�������ź�','�в��ź�')
title('��ȡ�źŵ�ѭ��Ƶ����ʾ');

% % �����źź���ȡ�źŵ�Ƶ�׷���
% Nw =  Fs;
% Nv =  ceil(3/4*Nw);
% Nfft = Nw;
% figure
% [S,f] = pwelch(xhat(:,end),Nw,Nv,Nfft,Fs);
% plot(f(1:5000),10*log10(S(1:5000))),hold on
% [S,f] = pwelch(y(:,end),Nw,Nv,Nfft,Fs);
% plot(f(1:5000),10*log10(S(1:5000)))
% % [S,f] = pwelch(nhat(:,1),Nw,Nv,Nfft,Fs);hold on
% % plot(f(1:5000),10*log10(S(1:5000)))
% xlabel('Ƶ�� (Hz)'),ylabel('dB')
% legend('��ȡ������','��������')
% title('�źŵ�Ƶ�׷���');

% figure
% [S,f] = pwelch(nhat(:,1),Nw,Nv,Nfft,Fs);hold on
% plot(f(1:5000),10*log10(S(1:5000)))
% title('�źŵ�Ƶ�׷���');