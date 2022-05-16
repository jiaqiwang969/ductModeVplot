%synthesized signal by noise signal + sum of sinusoidal signals
%STFT_LE_FD algorithm evaluation

clc
clear

A = 1;
fs = 1200;
N = 60000;
% adding pure sine components with closely spaced frequencies
periodic_sig = zeros(1,N);
n = [1:N]/fs;
f1 = 200;    %starting frequency
f2 = 5;   %frequency interval
Nf = 5;  %number of frequencies
for i_period = 1:1:Nf
    periodic_sig = periodic_sig + A*sin(2*pi*(f1+f2*i_period)*n);
end
%adding end

ratio = -5;  %信噪比
ratio = 10^(ratio/20);
amp = A/ratio;
randn('state',0);
noise = amp*randn(1,N);
data = periodic_sig + noise;
x = (data - mean(data))/std(data);
[row,col] = size(x);
if row<col
    x = x';
end

%*************************************************************
%STFT_LE
delay = 100;
Nwind = 256;
NFFT = 2*Nwind;
% version with absolute magnitude TF
G = STFT_LE(x,Nwind+delay,Nwind,round(2*Nwind/3),NFFT,1,'hanning');
G1 = abs(G);
[y1,g] = Filt_STFT(x,G1.');    %复数转置G1.'; 复数共轭转置G1';
N1 = length(y1);
e1 = x(1:N1,:) - y1;

%*************************************************************
%STFT_LE for compare
delay = 100;
Nwind = 2048;
NFFT_CMP = 2*Nwind;
% version with absolute magnitude TF
G_CMP = STFT_LE(x,Nwind+delay,Nwind,round(2*Nwind/3),NFFT_CMP,1,'hanning');
G1_CMP = abs(G_CMP);
[y1_CMP,g_CMP] = Filt_STFT(x,G1_CMP.');    %复数转置G1.'; 复数共轭转置G1';
N1 = length(y1_CMP);
e1_CMP = x(1:N1,:) - y1_CMP;

%

figure
Nfft = length(G);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
plot(f,G1(1:Nfft/2))
hold on
Nfft = length(G_CMP);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
plot(f,G1_CMP(1:Nfft/2),'r')
xlabel('G and G CMP')

figure,plot(n(1:N),x),hold on,plot(n(1:N),e1,'r'),plot(n(1:N),e1_CMP,'g')
xlabel('x, e and e CMP')
figure,plot(n(1:N),x),hold on,plot(n(1:N),y1,'r'),plot(n(1:N),y1_CMP,'g')
xlabel('x, y and y CMP')
[psd_x,f_x1] = pwelch(x,[],[],[],fs);
[psd_e1,f_e1] = pwelch(e1,[],[],[],fs);
[psd_e1_CMP,f_e1_CMP] = pwelch(e1_CMP,[],[],[],fs);
figure
plot(f_x1,10*log10(psd_x));
hold on
plot(f_e1,10*log10(psd_e1),'r');
plot(f_e1_CMP,10*log10(psd_e1_CMP),'g');
xlabel('x, e and e CMP')
[psd_y1,f_y1] = pwelch(y1,[],[],[],fs);
[psd_y1_CMP,f_y1_CMP] = pwelch(y1_CMP,[],[],[],fs);
figure
plot(f_x1,10*log10(psd_x));
hold on
plot(f_y1,10*log10(psd_y1),'r');
plot(f_y1_CMP,10*log10(psd_y1_CMP),'g');
xlabel('x, y and y CMP')



