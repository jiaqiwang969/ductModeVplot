function x = create_bearing_signal_noise(N,A,fs,fn,fp,var_tao,fr,snr)
% N         ���ݶ��ڷ��������������.(������ת�˶��ٸ����ڣ�ȡ������)
% A         �����ֵ: ���ݹ���λ�ò�ͬ, ����ܵ��ĵ�������.
%           ���� (ORF), �������ź� (IRF/REF).
% fs        ����Ƶ��.
% fn        ϵͳ����Ƶ��.
% fp        ��������Ƶ��.
% var_tao	����������ڹ����Ļ���.
% fr        תƵ. ���Ϊ0, ����ORF; ����ΪIRF/REF.

B = 800;
M = round(fs/fp);   % ������תһ������ʱ�Ĳ�������(tp/ts)
t = (0:M*N-1)*1/fs; % ʱ�䡣(M*N��������)
h = (0:M-1)*1/fs;   % ʱ�䣨������תһ�������ڵĲ����㣩

amp = A*(cos(2*pi*fr*t+2*pi*0.15)+1);
randn('state',0);%������������ӣ�������ͬ���������ͬ
tao = randn(1,N);
tao = round(tao*var_tao*M);
tao(1,1) = 0;
deltas = zeros(size(t));

for i = 1:N
    deltas(1+(i-1)*M+tao(1,i)) = amp(1+(i-1)*M+tao(1,i)); %ֻ�г������ʱ����������������壬LT
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