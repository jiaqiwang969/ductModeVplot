%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ѭ��ά���˲���+�������˲���
% �����������¶���������������ȡ��ͬʱ�����г���Ϳ������
% ��Ҫ�ο����ף�
% (1) Liang Yu, Haijun Wu, Jerome Antoni, Weikang Jiang��
% Extraction and imaging of aerodynamically generated sound field of rotor 
% blades in the wind tunnel test, Mechanical Systems and Signal Processing, 
% Volume 116, 1 February 2019, Pages 1017-1028. 

% (2) Liang Yu, Jerome Antoni, Haijun Wu, Weikang Jiang, 
% Reconstruction of cyclostationary sound source based on a back-propagating
% cyclic wiener filter, Journal of Sound and Vibration,2019,442 :787-799.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
clear
close all

SampeFreq = 51200 ;  % ����Ƶ��
% AnaChannelNums = 13;  % ������10ͨ��13 24
SampleTime = 10; % ����ʱ��
Fs = SampeFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ѭ��ά���˲���������������������ȡ������1�׺�2��ͳ����Ϣ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��ȡת���ļ�
% fid = fopen([DateFilePath,DateFileName,'-C-period.dat'], 'r');
% time_blade = fscanf(fid,'%f'); %��������ѹ�źţ�ÿһȦ��ʱ��
% fclose(fid);
if 1
    [Datafilename,mypath1]=uigetfile('D:\Aѧϰ����\����������\��������3\Bo105\ƽ��\*.*','���뽰Ҷ��λ����');
    S=strcat(mypath1,Datafilename);
    
    fid = fopen(S,'r','b');
    Nona=fread(fid,1,'single');
    Timenummic=fread(fid,1,'single');
    CLenthTime=fread(fid,1,'single');
    NLenthTime=fread(fid,1,'single');
    BladePositionSig1=fread(fid,CLenthTime*NLenthTime,'single');
    fclose(fid);
    
    
    Npiont=CLenthTime*NLenthTime/128;
%     Npiont=CLenthTime*NLenthTime
    for NNp=1:Npiont
        BladePositionSig(NNp)=sum(BladePositionSig1(((NNp-1)*128+1):NNp*128));
    end
    disp(['��ɷ�λ�����ݶ���!'])
    time_blade = BladePositionSig;
end 

% BO105 ���źŶ���
if 0
[Datafilename,mypath1]=uigetfile('D:\Aѧϰ����\����������\��������3\Bo105\ƽ��\*.*','���뽰Ҷ��λ����'); 
S=strcat(mypath1,Datafilename);

fid = fopen(S,'r','b');
Nona=fread(fid,1,'single');
Timenummic=fread(fid,1,'single');
CLenthTime=fread(fid,1,'single');
NLenthTime=fread(fid,1,'single');
BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
fclose(fid);
disp(['��ɷ�λ�����ݶ���!'])
time_blade = BladePositionSig;
end


% ��ȡ��������
[Datafilename,mypath1]=uigetfile('D:\Aѧϰ����\����������\��������3\Bo105\ƽ��\*.*','������������'); 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ��ԭʼ���ݼ�ʱƵ����
AnaChannelNums = 10;  % ������10ͨ��13 24
y = RawData(1:SampleTime*SampeFreq,AnaChannelNums); % ��ԭʼ�����ж�ȡ��������
L = size(y,1);  % �������ݵĳ���
Fs = SampeFreq;

% noise = y;
% signal = y = RawData(1:SampleTime*SampeFreq,1:10);
% SNR = 10*log10(var(signal)/var(noise))
% Lp = 20*log10(rms(signal)/(2*10^(-5)));

% ��ʾ�������ݺͻ�����Ƶ�׷���(ͼ1)
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
time_disp = 0.1;
frequency_disp = 6000;
% figure
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% set(findobj('type','axes'),'fontsize',28);
% xlabel('ʱ�� (��)'),ylabel('����'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('��������');

figure
subplot(211)
plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'k','LineWidth',0.5);
% xlabel('ʱ�� (��)'),ylabel('����'),
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);
% title('��������');
subplot(212)
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'k','LineWidth',1) % ֻ��ʾ��6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
% title(['���������ź�Ƶ�� ; Ƶ�ʷֱ��� = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');


% ʱƵ�׷���
Nw =  256;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
figure
spectrogram(y,Nw,Nv,Nfft, Fs,'yaxis')
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);

% ѭ��Ƶ������Է���
alpha_max = 300; % ����ѭ��Ƶ�� 300
Nw = 256; % ���ڳ��� 256
opt.coh = 1;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);

figure,
plot(alpha(2:end),sum(abs(S(:,2:end))),'LineWidth',2),
% xlabel('ѭ��Ƶ�� (Hz)','fontsize',12,'fontweight','b')
xlabel('Cyclic frequency (Hz)','fontsize',12,'fontweight','b')
ylabel('Sum of the spectral frequency (Hz)','fontsize',12,'fontweight','b');
% ylabel('��Ƶ����� (Hz)','fontsize',12,'fontweight','b');
% title('ѭ��Ƶ����ʾ');
set(findobj('type','axes'),'fontsize',12);
% xlim([0 alpha_max])

% % %  ʱ��ͬ��ƽ����ȡһ��ѭ��ƽ���ź�
% % % [~, syn_c_num] = size(time_blade);
% % Period_time = mean(time_blade);
% % Period_sample = round(Fs*Period_time*6);
% % Period_number = fix(L/Period_sample);
% % y_Period = y(1:Period_number*Period_sample);
% % y_Period_matrix = reshape(y_Period,Period_sample,Period_number);
% % My_Period = mean(y_Period_matrix,2);
% % ave_syn_y = repmat(My_Period,Period_number,1);
% % y_residual = y_Period - ave_syn_y;

% time_disp = 0.1;
% figure;
% subplot(311)
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('�����ź�')
% subplot(312)
% plot([0:1/Fs:time_disp-1/Fs],ave_syn_y(1:time_disp*Fs,end),'LineWidth',0.5);
% % xlabel('ʱ�� (��)'),ylabel('����'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('��ȡ��1��ѭ��ƽ���ź�')
% subplot(313)
% plot([0:1/Fs:time_disp-1/Fs],y_residual(1:time_disp*Fs,end),'LineWidth',0.5);
% % xlabel('ʱ�� (��)'),ylabel('����'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('�в��ź�')
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');

time_disp = 0.5;
figure
plot([0:1/Fs:time_disp-1/Fs],ave_syn_y(1:time_disp*Fs,end),'LineWidth',0.5);

Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
time_disp = 0.1;
frequency_disp = 5000;
figure
subplot(311)
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % ֻ��ʾ��6000Hz
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['�����ź�Ƶ�� ; Ƶ�ʷֱ��� = ',num2str(diff(f(1:2))),'Hz'])

subplot(312)
[S,f] = pwelch(ave_syn_y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % ֻ��ʾ��6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['��ȡ��1��ѭ��ƽ���ź�Ƶ�� ; Ƶ�ʷֱ��� = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');

subplot(313)
[S,f] = pwelch(y_residual,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % ֻ��ʾ��6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['�в��ź�Ƶ�� ; Ƶ�ʷֱ��� = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');

y = y_residual;
L = length(y);
%%  ѭ��ά���˲���
Na = 10; % ѭ��Ƶ�ʵĽ���
orders = (1:Na/2); % ѭ��Ƶ�ʵĸ���
orders = [orders;-orders]; 
orders = orders(:); % �ܵ�ѭ��Ƶ�ʸ���

alpha = 34.2/Fs; %34.2/Fs;%136.8/Fs/4; % ��һ��ѭ��Ƶ�� 136.8/Fs
phi = 2*pi*alpha*(0:L-1)'; % ��λ������
fi = alpha*ones(L,1); % Ƶ�ʵ�����

% ����ѭ�����׾���
Nw = 2^8;
[Syya,Syaya] = CyclicSpecMat(y,phi,fi,orders,Nw);
f = Syaya.f*Fs;

% ����ѭ��ά���˲���
Ns = 1; % ѭ��ƽ����Դ�ĸ���
[G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,1e-6);

% ����Gabor�任�����˲�
[xhat,cx] = CyclicSpatialFiltering(y,G(1,:,:),phi,fi,orders); %%xhat������Ŀ���ź�
[L_xhat,~] = size(xhat);
nhat = y(1:L_xhat) - xhat;  %%nhat��ʣ��������ź�

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ǰһ��Ĳ�����������ȡ�����ݽ��жԱ�
time_disp = 0.1;% 0.05;
figure;
% subplot(211);
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',2);
% % xlabel('ʱ�� [��])'),ylabel('����[Pa]');
% ylabel('Magnitude [Pa]','fontsize',12,'fontweight','b');
% xlabel('Time [Second]','fontsize',12,'fontweight','b');
% title('�в��ź�');
subplot(212);
plot([0:1/Fs:time_disp-1/Fs],xhat(1:time_disp*Fs,end),'LineWidth',2);
% xlabel('ʱ�� [��]'),ylabel('����[Pa]');
ylabel('Magnitude [Pa]','fontsize',12,'fontweight','b');
xlabel('Time [Second]','fontsize',12,'fontweight','b');
title('��һ����ȡ�Ķ���ѭ��ƽ���ź�');
% subplot(313);
% plot([0:1/Fs:time_disp-1/Fs],nhat(1:time_disp*Fs,end),'LineWidth',2);
% % xlabel('ʱ�� (��)'),ylabel('����');
% ylabel('Magnitude [Pa]','fontsize',12,'fontweight','b');
% xlabel('Time [Second]','fontsize',12,'fontweight','b');
% title('�в��ź�ʱ������');
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');

% % FFT psd
% frequency_disp = Fs/2;  %��ʾƵ�ʵ�**Hz
% f = Fs*(0:(L/2))/L;
% Y = fft(xhat); %FFT
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

% Nw =  Fs;   % window
% Nv =  ceil(3/4*Nw); % overlap
% Nfft = Nw;  % nfft
% frequency_disp = 6000;
% figure
% [S,f] = pwelch(residual,Nw,Nv,Nfft,Fs);
% plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'k','LineWidth',1.25) 
% ylabel('PSD [dB/Hz]','fontsize',12);
% xlabel('Frequency [Hz]','fontsize',12);
% title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
% xlim([0 frequency_disp])
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% set(gcf,'position',[400 250 500 250])

% figure
% plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% set(findobj('type','axes'),'fontsize',28);
% xlabel('ʱ�� (��)'),ylabel('����'),
% ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
% xlabel('Time [Second]','fontsize',12);
% title('��������');
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
time_disp = 0.1;
frequency_disp = 6000;
plot([0:1/Fs:time_disp-1/Fs],xhat(1:time_disp*Fs,end),'k','LineWidth',0.5);
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);
figure
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',1) % ֻ��ʾ��6000Hz
% xlabel('Frequency (Hz)'),ylabel('dB'),
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
% title(['���������ź�Ƶ�� ; Ƶ�ʷֱ��� = ',num2str(diff(f(1:2))),'Hz'])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

%% ���ݷ���
% �˲����ѭ��Ƶ������Է���
alpha_max = 300; % ����ѭ��Ƶ��
Nw = 2^8; % ���ڳ���
opt.coh = 1;

figure
[S,alpha,f,STFT,t,Nv] = Fast_SC(xhat(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end))),'k','LineWidth',2), hold on;
[S,alpha,f,STFT,t,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);
plot(alpha(2:end),sum(abs(S(:,2:end))),'r--','LineWidth',1);
% [S,alpha,f,STFT,t,Nv] = Fast_SC(nhat(:,end),Nw,alpha_max,Fs,opt);
% plot(alpha(2:end),sum(abs(S(:,2:end))),'y','LineWidth',2);
% xlabel('ѭ��Ƶ�� (Hz)')
xlabel('Cyclic frequency (Hz)','fontsize',12,'fontweight','b')
ylabel('Sum of the spectral fruequcy (Hz)','fontsize',12,'fontweight','b');
% legend('��ȡ���ź�','�������ź�','�в��ź�')
legend({'Extracted signal','Measured signal'},'FontSize',12)
% title('��ȡ�źŵ�ѭ��Ƶ����ʾ');
set(findobj('type','axes'),'fontsize',12);

% �����źź���ȡ�źŵ�Ƶ�׷���
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
frequency_disp = 800;
figure
[S,f] = pwelch(xhat(1:2*Fs),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',2),hold on
[S,f] = pwelch(y(1:2*Fs),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'r','LineWidth',2),hold on
xlabel('Ƶ�� (Hz)'),ylabel('dB')
ylabel('PSD (dB/Hz)','fontsize',12,'fontweight','b');
xlabel('Time [Second]','fontsize',12,'fontweight','b');
legend('��ȡ������','��������')
title('�źŵ�Ƶ�׷���');
set(findobj('type','axes'),'fontsize',12);

% ʱƵ�׷���
Nw =  256;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
figure
spectrogram(xhat,Nw,Nv,Nfft, Fs,'yaxis')
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  �������˲����������г���Ϳ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��ȡת���ļ�
% fid = fopen([DateFilePath,DateFileName,'-C-period.dat'], 'r');
% time_blade = fscanf(fid,'%f'); %��������ѹ�źţ�ÿһȦ��ʱ��
% fclose(fid);
if 1
    [Datafilename,mypath1]=uigetfile('E:\reseach\��������3\*.*','���뽰Ҷ��λ����');
    S=strcat(mypath1,Datafilename);
    
    fid = fopen(S,'r','b');
    Nona=fread(fid,1,'single');
    Timenummic=fread(fid,1,'single');
    CLenthTime=fread(fid,1,'single');
    NLenthTime=fread(fid,1,'single');
    BladePositionSig1=fread(fid,CLenthTime*NLenthTime,'single');
    fclose(fid);
    
    
    Npiont=CLenthTime*NLenthTime/128;
%     Npiont=CLenthTime*NLenthTime
    for NNp=1:Npiont
        BladePositionSig(NNp)=sum(BladePositionSig1(((NNp-1)*128+1):NNp*128));
    end
    disp(['��ɷ�λ�����ݶ���!'])
    time_blade = BladePositionSig;
end 

% BO105 ���źŶ���
if 0
[Datafilename,mypath1]=uigetfile('D:\��������3\*.*','���뽰Ҷ��λ����'); 
S=strcat(mypath1,Datafilename);

fid = fopen(S,'r','b');
Nona=fread(fid,1,'single');
Timenummic=fread(fid,1,'single');
CLenthTime=fread(fid,1,'single');
NLenthTime=fread(fid,1,'single');
BladePositionSig=fread(fid,CLenthTime*NLenthTime,'single');
fclose(fid);
disp(['��ɷ�λ�����ݶ���!'])
time_blade = BladePositionSig;
end

% ����ת����Ϣͬ��һ��˲ʱƵ��
f_rot=[]; % ˲ʱƵ��
for i=1:length(time_blade)
f_i=1/time_blade(i); % ��Ƶ��1/time_blade(i) ��ҶƵ��4�ǽ�Ҷ�ĸ���
f_rot=[f_rot f_i*ones(1,round(time_blade(i)*Fs))];
end
f_rot = f_rot.';

% �ź�׼��
ratio = 2;
y_xhat = xhat(1:ratio*Fs); 
% y_xhat = xhat(1:2*Fs);   %ȡ�����ź��ͽ��������˲��������˲�
fp = f_rot(1:size(y_xhat,1)); 
% t = (1:size(y_xhat,2))/Fs;                                 %��ɢʱ���ź�

%% �������˲�������׼��
r = 2000;   % ����Ȩ�� 2000
ford = 1;    % �����˲�������Ϊ1��Ҳ����ȡ2��
tol = 0.1;    % ����
maxit = 500;  % ����������

V_fp = [4*fp 2*4*fp 3*4*fp 4*4*fp 5*4*fp 6*4*fp]; % 4 ��ҶƬ�ĸ���, 5 ��ϣ���������г������
[nt,nch] = size(V_fp);

[xm,bwm,Tm,fl,rr,it,rv,xr] = vkm(y_xhat,V_fp,Fs,r*ones(nch,1),ford,tol,maxit);
figure,plot(0:maxit, rv,'LineWidth',2);
ylabel('residual error','fontsize',12,'fontweight','b');
xlabel('Iteration','fontsize',12,'fontweight','b');
%% ��ȡ��г���źŵı��
figure;
subplot(511)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,1),'LineWidth',2); 
legend({'ҶƵ��1г��'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(512)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,2),'LineWidth',2)
legend({'ҶƵ��2г��'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(513)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,3),'LineWidth',2)
legend({'ҶƵ��3г��'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(514)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,4),'LineWidth',2)
legend({'ҶƵ��4г��'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(515)
plot(0:1/Fs:0.1-1/Fs,xr(1:0.1*Fs,5),'LineWidth',2)
legend({'ҶƵ��5г��'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');

%% �����������лָ�г���źţ�Ƶ�ף�
Nw =  Fs;
Nv =  ceil(3/4*Nw);
Nfft = Nw;
frequency_disp = 1500;
figure
[S,f] = pwelch(y_xhat,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,1),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,2),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,3),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,4),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,6),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
hold on
[S,f] = pwelch(xr(:,5),Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'--','LineWidth',2)
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');
legend({'��Ƶ+г��','ҶƵ��1г��','ҶƵ��2г��','ҶƵ��3г��','ҶƵ��4г��','ҶƵ��5г��','ҶƵ��6г��'},'FontSize',12)
%%

% �����������лָ�г���źţ�Ƶ�ף�
% frequency_disp = 800;
% figure
% Spec_y = fft(y(1:ratio*Fs,1));
% f = linspace(0,Fs,length(Spec_y));
% plot(f,10*log10(abs(Spec_y/size(Spec_y,2))));
% xlim([0 frequency_disp])
% hold on
% Spec_xr1 = fft(xr(:,1));
% f = linspace(0,Fs,length(Spec_y));
% plot(f,10*log10(abs(Spec_xr1/size(Spec_xr1,2))));
% xlim([0 frequency_disp])
% Spec_xr2 = fft(xr(:,2));
% plot(f,10*log10(abs(Spec_xr2/size(Spec_xr2,2))));
% xlim([0 frequency_disp])
% Spec_xr3 = fft(xr(:,3));
% plot(f,10*log10(abs(Spec_xr3/size(Spec_xr3,2))));
% xlim([0 frequency_disp])
% Spec_xr4 = fft(xr(:,4));
% plot(f,10*log10(abs(Spec_xr4/size(Spec_xr4,2))));
% xlim([0 frequency_disp])
% Spec_xr5 = fft(xr(:,5));
% plot(f,10*log10(abs(Spec_xr5/size(Spec_xr5,2))));
% xlim([0 frequency_disp])


% ���������г���źźͿ���ź�
residual = y_xhat - sum(xr,2);
harmocs  = sum(xr,2);
figure,
subplot(211)
plot(0:1/Fs:0.1-1/Fs,harmocs(1:0.1*Fs))
legend({'ҶƵг������'},'FontSize',12);
ylabel('Magnitude','fontsize',12,'fontweight','b');
xlabel('time [second]','fontsize',12,'fontweight','b');
subplot(212)
plot(0:1/Fs:0.1-1/Fs,xhat(1:0.1*Fs),'LineWidth',0.5);
% legend({'��Ƶ����'},'FontSize',12);
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);

plot([0:1/Fs:time_disp-1/Fs],y(1:time_disp*Fs,end),'LineWidth',0.5);
% xlabel('ʱ�� (��)'),ylabel('����'),
ylabel('Magnitude [Pa]','fontsize',12); % ,'fontweight','b'
xlabel('Time [Second]','fontsize',12);
% % (3) �����׷�����Welch��s power spectral density estimate��
% Nw =  Fs/2;   % window
% Nv =  ceil(3/4*Nw); % overlap
% Nfft = Nw;  % nfft
% frequency_disp = Fs/4;
% figure
% [S,f] = pwelch(residual,Nw,Nv,Nfft,Fs);
% plot(f(1:frequency_disp/2),10*log10(S(1:frequency_disp2/)),'LineWidth',1) % ֻ��ʾ��6000Hz
% % xlabel('Frequency (Hz)'),ylabel('dB'),
% ylabel('PSD [dB/Hz]','fontsize',12);
% xlabel('Frequency [Hz]','fontsize',12);
% % title(['��ȡ��1��ѭ��ƽ���ź�Ƶ�� ; Ƶ�ʷֱ��� = ',num2str(diff(f(1:2))),'Hz'])
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% (3) �����׷�����Welch��s power spectral density estimate��
Nw =  Fs/2;   % window
Nv =  ceil(3/4*Nw); % overlap
Nfft = Nw;  % nfft
frequency_disp = Fs/8;
figure
[S,f] = pwelch(xhat,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp/2),10*log10(S(1:frequency_disp/2)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
% title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
ylim([-80 0])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 200])

% �����źš�г���źźͿ���źŵ�Ƶ��
frequency_disp = 2000;
figure
[S,f] = pwelch(y_xhat,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'b','LineWidth',2)
hold on
[S,f] = pwelch(harmocs,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'r--','LineWidth',2)
hold on
[S,f] = pwelch(residual,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'y--','LineWidth',2)
legend({'�����ź�','ҶƵг������','��Ƶ����'},'FontSize',12)
set(findobj('type','axes'),'fontsize',12);
ylabel('PSD [dB/Hz]','fontsize',12,'fontweight','b');
xlabel('Frequency [Hz]','fontsize',12,'fontweight','b');

