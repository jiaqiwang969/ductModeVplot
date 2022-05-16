clc;
clear;
close all;
N=30;
A=1;
fs=25600;
fn=3000;
fp=43;
var_tao=0.03;
Q=12;
% alpha=43;
M=8; 
L=1;
fr=0;
snr=-12;
x = create_bearing_signal_noise(N,A,fs,fn,fp,var_tao,fr,snr);
sizex=length(x);
figure;
subplot(3,1,1);
t=1/fs*(0:(sizex-1));
plot(t,x);
y1=2*abs(fft(x,sizex));
subplot(3,1,2);
f=fs/sizex*(0:(sizex-1));
plot(f,y1);
% plot(fs/sizex*(1:sizex/2.56),y1(1:sizex/2.56));
y2=2*abs(fftshift(fft(x,sizex)));
subplot(3,1,3);
plot(y2);
for i=0:1000
    alpha=i;
    [S,df]=cyc_spectral_crocorr_certainalpha(x,alpha,fs,'triang',M);
    if i==0
        S0=sum(abs(S.^2));
    end
    A(i+1)=sum(abs(S.^2))/S0;
end
figure;
plot(A)
