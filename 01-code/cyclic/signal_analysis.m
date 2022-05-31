%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%�źŷ���%%%%%%%%%%%%%%%
% ��Ҫ����������
% 1.������������
% (1) �ź�ʱ����
% (2) FFT Ƶ�׷���
% (3) PSD �����׷���
% (4) �����׷���
% (5) STFT ʱƵ����
% (6) Spectral Kurtosis/Kurtogram ���Ͷ�
% 2.ѭ��ƽ�ȷ���
% (7) SC ѭ�������/ѭ������ɷ���(Fast-SC�㷨)
% 3.ģ̬����
% (8) EMD ����ģ̬�ֽ�
%%%%%%%%%%%%��ӭ��������%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear,close all;
SampeFreq =51200 ;  % ����Ƶ��(Hz)
% SampleTime = 30; % ����ʱ��(s)
Fs = SampeFreq;
%% ��ȡ��������
[DataFileName,DataFilePath]=uigetfile('','ѡ��ԭʼ����'); % �򿪶�ȡ���ݴ���
load(strcat(DataFilePath,DataFileName));
disp('���ԭʼ���ݶ���!')
RawData =xhat;    % xΪ���������еı��������������RawData��
[m,n] = size(RawData);
if m<n
    RawData = RawData'; % ��RawDataͳһ��дΪ (��������*ͨ����) ����ʽ
end
AnaChannelNums = 30;  % ������ͨ����
AnaTime = 1;   % ������ʱ�䳤��(s)
y = RawData(1:AnaTime*SampeFreq,AnaChannelNums); % ��ԭʼ�����ж�ȡ��������
L = size(y,1);  % �������ݵĳ���
%% ������������
% (1) �������ݵ�ʱ����
% figure
% t = 0:1/Fs:AnaTime-1/Fs;    % ����ʱ�����
% plot(t,y,'k','LineWidth',1.5);
% time_disp = 10;    % ��ʾʱ�䳤��
% xlim([0 time_disp])
% ylabel('Amplitude');    % ע�⣺����ϵ�λ[m/s^2]��[Pa]
% xlabel('Time [Second]');
% title('Time Domain Waveform','fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(findobj('type','axes'),'fontweight','b');
% set(gcf,'position',[400 250 500 250])
% % 
% % (2) FFTƵ�׷���
frequency_disp = Fs/2;  %��ʾƵ�ʵ�**Hz
f = Fs*(0:(L/2))/L;
Y = fft(y); %FFT
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure
plot(f,P1,'k','LineWidth',1.5) 
title('FFT Amplitude Spectrum','fontweight','b')
xlabel('f [Hz]')
ylabel('Amplitude')    % ע�⣺����ϵ�λ[m/s^2]��[Pa]
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])


% (3) �����׷�����Welch��s power spectral density estimate��
Nw =  Fs;   % window
Nv =  ceil(3/4*Nw); % overlap
Nfft = Nw;  % nfft
frequency_disp = Fs/2;
figure
[S,f] = pwelch(y,Nw,Nv,Nfft,Fs);
plot(f(1:frequency_disp),10*log10(S(1:frequency_disp)),'k','LineWidth',1.25) 
ylabel('PSD [dB/Hz]','fontsize',12);
xlabel('Frequency [Hz]','fontsize',12);
title(['Signal PSD (Frequency resolution=',num2str(diff(f(1:2))),'Hz)'],'fontweight','b')
xlim([0 frequency_disp])
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

% (4) �����׷���
frequency_disp = 2e3;  %��ʾƵ�ʵ�2000Hz
f = Fs*(0:(L/2))/L;
y1=y-mean(y);   %ȥ��ֱ������
Yh=abs(hilbert(y1));    %ϣ�����ر任
lnY=length(Yh);
Y1=fft(Yh-mean(Yh))/lnY;    %�����FFT
Y1=abs(Y1(1:L/2+1));
figure
plot(f,Y1,'k','LineWidth',1.5);
title('Envelop Spectrum','fontweight','b')
xlim([0 frequency_disp]) %Ƶ�ʵ���ʾ��Χ
xlabel('f [Hz]')
ylabel('Amplitude')    % ע�⣺����ϵ�λ[m/s^2]��[Pa]
set(findobj('type','axes'),'fontsize',12);
set(findobj('type','axes'),'fontweight','b');
set(gcf,'position',[400 250 500 250])

% 
% % (5)STFTʱƵ�׷�����Spectrogram using short-time Fourier transform��
frequency_disp = Fs/2;    % ��ʾƵ�ʳ���
Nw =  2^11;  % ����
window = hanning(Nw);      % hanning��
Nv =  ceil(3/4*Nw); % overlap:�ص�������:3/4*Nw
Nfft = 4*Nw;  % nfft:DFT����
[S, F, T] = spectrogram(y, window, Nv, Nfft, Fs);
figure
imagesc(T, F, log10(abs(S)))
set(gca, 'YDir', 'normal')
xlabel('Time [Seconds]')
ylabel('Frequency [Hz]')
title('STFT Spectrum','fontweight','b')
h=colorbar;
set(get(h,'Title'),'string','dB/Hz');
ylim([0 frequency_disp]) %Ƶ�ʵ���ʾ��Χ
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);
set(gcf,'position',[400 250 500 350])
% 
% % (6) Spectral Kurtosis/Kurtogram ���Ͷ�
% figure
% pkurtosis(y,Fs) % Spectral Kurtosis
% title('Spectral Kurtosis','fontweight','b')
% ylabel('Spectral Kurtosis'),xlabel('f [kHz]'),
% xlim([0 Fs/2000])
% set(findobj('type','axes'),'fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(gcf,'position',[400 250 500 250])
% 
% figure
% kurtogram(y,Fs) % Kurtogram
% xlim([0 Fs/2000])
% set(findobj('type','axes'),'fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(gcf,'position',[400 250 500 400])
%% ѭ��ƽ�ȷ���
% (7) SC ѭ�������/ѭ������ɷ���  (Fast-SC�㷨)
% ����������Ӻ��� Fast_SC.m
alpha_max = 800; % ����ѭ��Ƶ��
Nw = 2^8; % ���ڳ���
opt.coh = 1;    %Ϊ0ʱ����Spectral Correlation��Ϊ1ʱ����Spectral Coherence
[S,alpha,f,STFT,~,Nv] = Fast_SC(y(:,end),Nw,alpha_max,Fs,opt);  %����Fast-SC�㷨����SC

% ���������/�����
figure
surf(alpha(2:end),f,abs(S(:,2:end))),axis xy,colorbar()
shading flat
if opt.coh == 0,title('Spectral Correlation','fontweight','b'),else,title('Spectral Coherence','fontweight','b'),end
xlabel('cyclic frequency \alpha [Hz]'),ylabel('f [Hz]'),
xlim([0 alpha_max]),ylim([0 Fs/2])
view(30,30)
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);
set(gcf,'position',[400 250 500 350])
% ����ƽ�������/ƽ�������
figure
plot(alpha(2:end),mean(abs(S(:,2:end))),'k'),xlim([0 alpha_max])
if opt.coh == 0,title('Mean Spectral Correlation','fontweight','b'),else,title('Mean Spectral Coherence','fontweight','b'),end
xlabel('cyclic frequency \alpha (Hz)')
set(findobj('type','axes'),'fontweight','b');
set(findobj('type','axes'),'fontsize',12);
set(gcf,'position',[400 250 500 250])
% %% ģ̬����
% % (8) EMD ����ģ̬�ֽ�
% [imf_emd,residual]=emd(y,'MaxNumIMF',5);
% imf_emd = imf_emd';
% figure;
% frequency_disp = Fs/2;  %��ʾƵ�ʵ�**Hz
% subplot(size(imf_emd,1)+1,2,1);
% t = 0:1/Fs:AnaTime-1/Fs;    % ����ʱ�����
% plot(t,y,'k','LineWidth',1.5);
% title('EMD');
% ylabel('Signal')
% subplot(size(imf_emd,1)+2,2,2);
% f = Fs*(0:(L/2))/L;
% Y = fft(y); %FFT
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% plot(f,P1,'k','LineWidth',1.5);
% xlim([0 frequency_disp]) %Ƶ�ʵ���ʾ��Χ
% title('Amplitude Spectrum');
% ylabel('Signal')
% for i = 2:size(imf_emd,1)+1
%     subplot(size(imf_emd,1)+2,2,i*2-1);
%     plot(t,imf_emd(i-1,:),'k','LineWidth',1.5);
%     ylabel(['IMF',num2str(i-1)])
%     subplot(size(imf_emd,1)+2,2,i*2);
%     Y = fft(imf_emd(i-1,:)); %FFT
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     plot(f,P1,'k','LineWidth',1.5);
%     ylabel(['IMF',num2str(i-1)])
%     xlim([0 frequency_disp]) %Ƶ�ʵ���ʾ��Χ
% end
% subplot(size(imf_emd,1)+2,2,2*size(imf_emd,1)+3);
% plot(t,residual,'k','LineWidth',1.5);
% xlabel('Time [s]')
% ylabel('Residual')
% subplot(size(imf_emd,1)+2,2,2*size(imf_emd,1)+4);
% Y = fft(residual); %FFT
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% plot(f,P1,'k','LineWidth',1.5);
% xlabel('Frequency [Hz]')
% ylabel('Residual')
% xlim([0 frequency_disp]) %Ƶ�ʵ���ʾ��Χ
% set(findobj('type','axes'),'fontweight','b');
% set(findobj('type','axes'),'fontsize',12);
% set(gcf,'position',[300 -400 800 850])