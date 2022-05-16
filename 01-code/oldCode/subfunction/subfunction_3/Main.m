clc
clear
N=5120;
f=5120;
n=0:(N-1);
t=n/f;
x=sin(2*pi*50*t);
X=fft(x);
f=f*(0:N-1)/N;     % 被修改
subplot(2,1,1);
plot(f,abs(X)*2/N,'-');     % 被修改
axis([0,115,0,1.2]);     % 被修改
set(gca, 'XTickMode', 'manual', 'XTick', [0, 20,30,40,45,50,55,60,70,80,90,100,110,115]);
xlabel('Frequency');
ylabel('|F(k)|');
grid on;

f0=45;%移动的频率
h=exp(-j*2*pi*f0*t);
y=x.*h;%数字变频，把w=100点移动到原点
y1=resample(y,1,10);%重新采样，采样频率为fs/N
k1=0:1:511;
h1=exp((j*2*pi*f0*k1)/512);
y1=y1.*h1;%反移频，使得前后频率一致
y1=abs(fft(y1));%快速傅里叶变换
subplot(2,1,2);
plot(k1,y1*2/512);     % 被修改
axis([45,55,0,1.2]);     % 被修改
set(gca, 'XTickMode', 'manual', 'XTick', [45,46,47,48,49,50,51,52,53,54,55]);
title('细化10倍后的频率特性');
xlabel('频率');
ylabel('幅度值');
grid on;
A=features_time_domain(x)