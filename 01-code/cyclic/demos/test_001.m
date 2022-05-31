clear;
%
addpath '..';
%
fs=120;%采样频率
nt=1200;%采样点数
t=(1:nt)./fs;%时间序列
t=t(:);
%信号1频率
f0=20;f1=40;
opm=log2(f1/f0)*60/t(end);
fup=f0*2.^(opm/60*t);fup=fup(:);%频率从f0-f1
%信号2频率
fc=30*ones(size(fup));%频率为fc
%幅值
am1=0.5+0.5/t(end)*t;
am2=1*ones(size(am1));
%构建两个正弦信号
x1=am1.*sin(2*pi*cumsum(fup/fs));
x2=am2.*sin(2*pi*cumsum(fc/fs));
%
fp=[fc,fup];
xp=[x2,x1];
[nt,nch]=size(xp);
%
nsens=1;
%A=randn(nsens,nch);
A=[1,1];
%
xe=(A*xp.').';
xc=A(ones(nt,1),:).*xp;
%
%xe=xe+0.05*randn(nt,nsens); % add noise
%
nfft=128;nwin=100;novlap=50;dflag='mean';
[Sxx_m,freq,time]=p_gram(xe,nfft,fs,nwin,novlap,dflag);
%
r=1600;
ford=1;
[xs(:,1),bw,T,xr(:,1)] = vk2(xe,fp(:,1),fs,r,ford);
[xs(:,2),bw,T,xr(:,2)] = vk2(xe,fp(:,2),fs,r,ford);
%
figure(1);
h=plot(t,xc,'--',t,abs(xs),'r');
set(h(2:end),'linewidth',2);
xlabel('Time [Sec]');
ylabel('Amplitude');
title('Single Order Extraction');
%
rm=1600;
tol=0.001;maxit=1000;
[xm,bwm,Tm,fl,rr,it,rv,xrm] = vkm(xe(:,1),[fc,fup],fs,rm*ones(nch,1),ford,tol,maxit);
%
figure(2);
h=plot(t,xc,'--',t,abs(xm),'r');
set(h(2:end),'linewidth',2);
xlabel('Time [Sec]');
ylabel('Amplitude');
title('Multiple Order Extraction');
%
figure(3);
semilogy(rv);
xlabel('Iteration');
ylabel('Residual');
%
figure(4);
typ='lin';f_min=0;f_max=50;
spec_plot2(Sxx_m,freq,time,typ,f_min,f_max,0,0.2);
colorbar('off');
g=colormap('gray');colormap(flipud(g));
title('Spectrogram of orders');
%
figure(5);
f = fs*(0:(nt/2))/nt;   % Define the frequency domain f
Y_xe = fft(xe);
P2_xe = abs(Y_xe/nt);
P1_xe = P2_xe(1:nt/2+1);
P1_xe(2:end-1) = 2*P1_xe(2:end-1);
subplot(311);
plot(f,P1_xe,'k');
legend('Original Signal');
xlabel('f (Hz)')
ylabel('Amplitude')
title('Amplitude Spectrum of Original Signal')
%
Y_xr1 = fft(xr(:,1));
P2_xr1 = abs(Y_xr1/nt);
P1_xr1 = P2_xr1(1:nt/2+1);
P1_xr1(2:end-1) = 2*P1_xr1(2:end-1);
Y_xr2 = fft(xr(:,2));
P2_xr2 = abs(Y_xr2/nt);
P1_xr2 = P2_xr2(1:nt/2+1);
P1_xr2(2:end-1) = 2*P1_xr2(2:end-1);
subplot(312);
plot(f,P1_xr1,'r',f,P1_xr2,'b');
legend('Extraction1','Extraction2');
title('Amplitude Spectrum of Single Order Extraction')
xlabel('f (Hz)')
ylabel('Amplitude')
%
Y_xrm1 = fft(xrm(:,1));
P2_xrm1 = abs(Y_xrm1/nt);
P1_xrm1 = P2_xrm1(1:nt/2+1);
P1_xrm1(2:end-1) = 2*P1_xrm1(2:end-1);
Y_xrm2 = fft(xrm(:,2));
P2_xrm2 = abs(Y_xrm2/nt);
P1_xrm2 = P2_xrm2(1:nt/2+1);
P1_xrm2(2:end-1) = 2*P1_xrm2(2:end-1);
subplot(313);
plot(f,P1_xrm1,'r',f,P1_xrm2,'b');
legend('Extraction1','Extraction2');
title('Amplitude Spectrum of Multiple Order Extraction')
xlabel('f (Hz)')
ylabel('Amplitude')

figure(6)
subplot(211)
plot(t,x1,'b',t,x2,'r',t,xr(:,2),'m',t,xr(:,1),'k')
legend('x1','x2','Extraction1','Extraction2');
title('Single Order Extraction')
xlabel('t (s)')
ylabel('Amplitude')
subplot(212)
plot(t,x1,'b',t,x2,'r',t,xrm(:,2),'m',t,xrm(:,1),'k')
legend('x1','x2','Extraction1','Extraction2');
title('Multiple Order Extraction')
xlabel('t (s)')
ylabel('Amplitude')
figure(7)

plot(t,xe)
xlabel('t (s)')
ylabel('Amplitude')
title('s(t)')
