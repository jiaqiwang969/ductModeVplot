%evaluation of theoretical frequency gain of ALE filter
%evaluation1 window type comparison: rectwin hanning parzenwin
clc
clear

Nwind = 512;
Nfft = 16*Nwind;
SNR = 1;

%rectwin
windowname = 'rectwin';
w = feval(windowname,Nwind);
w = w/sum(w);
W = fftshift(fft(w,Nfft));
H_abs = (0.5*SNR*Nwind*W.*conj(W))./(1+0.5*SNR*Nwind*W.*conj(W));
f = [-(Nfft/2):1:(Nfft/2-1)]/Nfft;
figure
plot(f,H_abs)
hold on
plot(f,abs(W),'r')
title('Comparison of frequency gain between ALE filter and its window')
xlabel(['window type: ',windowname])
axis([-0.02 0.02 0 1])

%hanning
windowname = 'hanning';
w = feval(windowname,Nwind);
w = w/sum(w);
W = fftshift(fft(w,Nfft));
H_abs = (0.5*SNR*Nwind*W.*conj(W))./(1+0.5*SNR*Nwind*W.*conj(W));
f = [-(Nfft/2):1:(Nfft/2-1)]/Nfft;
figure
plot(f,H_abs)
hold on
plot(f,abs(W),'r')
title('Comparison of frequency gain between ALE filter and its window')
xlabel(['window type: ',windowname])
axis([-0.02 0.02 0 1])

%parzenwin
windowname = 'parzenwin';
w = feval(windowname,Nwind);
w = w/sum(w);
W = fftshift(fft(w,Nfft));
H_abs = (0.5*SNR*Nwind*W.*conj(W))./(1+0.5*SNR*Nwind*W.*conj(W));
f = [-(Nfft/2):1:(Nfft/2-1)]/Nfft;
figure
plot(f,H_abs)
hold on
plot(f,abs(W),'r')
title('Comparison of frequency gain between ALE filter and its window')
xlabel(['window type: ',windowname])
axis([-0.02 0.02 0 1])

%diric -sanc
f = [-(Nfft/2):1:(Nfft/2-1)]/Nfft;
Diric_abs = abs(diric(2*pi*f,Nwind));
H_abs = (0.5*SNR*Nwind*Diric_abs)/(1+0.5*SNR*Nwind);
figure
plot(f,(H_abs),'k')
axis([-0.02 0.02 0 1])