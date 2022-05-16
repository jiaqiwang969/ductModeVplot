function [x,RSN] = create_bearing_signal(N,A,fs,fn,fp,var_tao,fr)

% N:        ���ݶ��ڷ��������������
% A:        �����ֵ: ���ݹ���λ�ò�ͬ������ܵ��ĵ������� ����Ϊ��������Ȧ���ϣ����������źţ���Ȧ���������ϣ�
% fs:       �źŵĲ���Ƶ�� 
% fn:       ϵͳ����Ƶ��
% fp:       ��������Ƶ��
% var_tao:  ��������ڹ������ڵĻ�����
% fr:       תƵ

B = 2000;
M=round(fs/fp);
t=(0:M*N-1)/fs;
h=(0:M-1)/fs;
                      

amp=A*(cos(2*pi*fr*t+2*pi*0.15)+1);
% amp=A*(cos(2*pi*fr*t+2*pi*0.15)+1+1*randn(1,length(t)));

deltas=zeros(size(t));    
randn('state',0);
tao=randn(1,N);
tao=round(tao*var_tao*M);             

tao(1,1)=0;
for i=1:N
    deltas(1+(i-1)*M+tao(1,i))=amp(1+(i-1)*M+tao(1,i));
end

s=exp(-1*B*h).*cos(2*pi*fn*h); 
tmpx=conv(deltas,s);
x_s=tmpx(1:M*N);


avr=1;
x=x_s;

if nargout==2
    RSN=20*log10(std(x_s)/avr);
elseif nargout==0
    figure(gcf+1);
    plot(t,x','k');
    title('signal');
    xlabel('time[s]');
    ylabel('amplitude');
end
