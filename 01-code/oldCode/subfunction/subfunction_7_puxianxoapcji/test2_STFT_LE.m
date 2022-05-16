% Test 2 experimental data from randall :

%  
clc
clear

d = 7.12;
D = 38.5;
n = 12;
fr = 6;
fo = (1-d/D)*n*fr/2;
fi = (1+d/D)*n*fr/2;

load fault3.mat;
data_fault = fault3(:,1);
data_fault = (data_fault - mean(data_fault))/std(data_fault);
N = length(data_fault);
fs = 48000;
load goodbear.mat;
data_good = goodbear(:,1);
data_good = (data_good - mean(data_good))/std(data_good);

[psd_good,f] = pwelch(data_good,[],[],1024,fs);
figure
plot(f,10*log10(psd_good));
hold on
[psd_fault,f] = pwelch(data_fault,[],[],1024,fs);
plot(f,10*log10(psd_fault),'r');

%band-pass filtration

f=[2*15000/fs 2*16000/fs 2*17000/fs 2*18000/fs];
%f=[2*15000/fs 2*16000/fs 2*17000/fs 2*18000/fs];(good for compare ALE and band envelope of fault 3)
A=[0 1 0];
rp=0.153;
rs=16.92;
devp=1-10^(-rp/20);
devs=10^(-rs/20);
dev=[devp devs devp];
[n,f0,A0,w]=remezord(f,A,dev);

if rem(n,2)
   n=n+1;
end

b=remez(n,f0,A0,w);

[h,w]=freqz(b,1,length(b),1);
hr=abs(h);
h=abs(h);
h=20*log10(h);
figure
subplot(211)
plot(b);grid;
subplot(212)
plot(w,h);grid;

% data_fault = data_good;
data_fault_filter = filter(b,1,data_fault);  
N = length(data_fault_filter);

%STFT_LE

delay = 100;
Nwind = 8192;
NFFT = 2*Nwind;

% version with absolute magnitude TF
G = STFT_LE(data_fault_filter,Nwind+delay,Nwind,round(2*Nwind/3),NFFT,1,'parzenwin');
G1 = abs(G);
[y1,g] = Filt_STFT(data_fault_filter,G1.');    %复数转置G1.'; 复数共轭转置G1';
N1 = length(y1);
e1 = data_fault_filter(1:N1,:) - y1;
% version with complex TF
G2 = G;
y2 = Filt_STFT(data_fault_filter,G2.');
N2 = Nwind + delay;
e2 = data_fault_filter(N2+1:N) - y2(1:N-N2);
    
t = (1:N)/fs;
figure,plot(t(1:N),data_fault_filter),hold on,plot(t(1:N),y1,'r'),plot(t(N2+1:N),y2(1:N-N2),':g')
figure,plot(t(1:N),data_fault_filter),hold on,plot(t(1:N),e1,'r'),plot(t(N2+1:N),e2(1:N-N2),':g')
[psd_x,f_x1] = pwelch(data_fault_filter);
[psd_e1,f_e1] = pwelch(e1);
[psd_e2,f_e2] = pwelch(e2);
figure
plot(f_x1,10*log10(psd_x));
hold on
plot(f_e1,10*log10(psd_e1),'r');
plot(f_e2,10*log10(psd_e2),':g');
[psd_y1,f_y1] = pwelch(y1);
[psd_y2,f_y2] = pwelch(y2);
figure
plot(f_x1,10*log10(psd_x));
hold on
plot(f_y1,10*log10(psd_y1),'r');
plot(f_y2,10*log10(psd_y2),':g');

% Demodulation

Nfft = length(data_fault_filter);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_data_fault_filter = hilbert(data_fault_filter);
envelope_spec_filter = abs(fft(abs(hilbert_data_fault_filter)));
figure
plot(f,envelope_spec_filter(1:Nfft/2)./Nfft)
hold on
Nfft = length(e1);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_e1 = hilbert(e1);
envelope_spec_e1 = abs(fft(abs(hilbert_e1)));
plot(f,envelope_spec_e1(1:Nfft/2)./Nfft,'r')
hold on
Nfft = length(e2);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_e2 = hilbert(e2);
envelope_spec_e2 = abs(fft(abs(hilbert_e2)));
plot(f,envelope_spec_e2(1:Nfft/2)./Nfft,'g')



