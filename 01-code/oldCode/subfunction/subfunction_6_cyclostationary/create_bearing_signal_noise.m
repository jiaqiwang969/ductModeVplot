function x = create_bearing_signal_noise(N,A,fs,fn,fp,var_tao,fr,snr)
% N         数据段内发生冲击的周期数.(故障旋转了多少个周期，取样长度)
% A         冲击幅值: 根据故障位置不同, 冲击受到的调制作用.
%           常数 (ORF), 或周期信号 (IRF/REF).
% fs        采样频率.
% fn        系统共振频率.
% fp        故障特征频率.
% var_tao	滚动体相对于滚道的滑动.
% fr        转频. 如果为0, 则是ORF; 否则为IRF/REF.

B = 800;
M = round(fs/fp);   % 故障旋转一个周期时的采样点数(tp/ts)
t = (0:M*N-1)*1/fs; % 时间。(M*N个采样点)
h = (0:M-1)*1/fs;   % 时间（故障旋转一个周期内的采样点）

amp = A*(cos(2*pi*fr*t+2*pi*0.15)+1);
randn('state',0);%生成随机数种子，种子相同，随机数相同
tao = randn(1,N);
tao = round(tao*var_tao*M);
tao(1,1) = 0;
deltas = zeros(size(t));

for i = 1:N
    deltas(1+(i-1)*M+tao(1,i)) = amp(1+(i-1)*M+tao(1,i)); %只有冲击发生时的随机滑动才有意义，LT
end
s = exp(-1*B*h).*cos(2*pi*fn*h);
x = conv(deltas,s);
x = x(1:M*N);
avr=std(x)/10^(snr/20);
noise=avr*randn(1,M*N);
x=x+noise;
if nargout == 0
    figure;
    plot(t,x,'k');
    title('signal');
    xlabel('time (s)');
    ylabel('amplitude');
end