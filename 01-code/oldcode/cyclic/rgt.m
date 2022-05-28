function [c,E] = rgt(x,g,R)
% c = rgt(x,g,R)
% ���ź�x����Gabor �任
% ���룺 
% x ��������źţ�
% gΪ�������ڳ��ȣ�
% RΪ���ں���ÿ���ƶ��ĳ��ȣ�
%����� 
% c ΪGabor �任���ϵ�����������е���Ŷ�ӦƵ�ʣ��е���Ŷ�Ӧʱ��
% E ΪGabor �任���ϵ���������ʱ����Ҷ�任���ϵ������֮�����λ��ϵ��c_Gabor = c_STFT.*E


x = x(:);
g = g(:);
L = length(x);
Nw = length(g);
if rem(Nw,2 > 0)
    error('the length of window g must be even !');
end
K = floor((L-Nw)/R+1);

c = zeros(Nw/2+1,K);
ind = 1;
wi = 2i*pi*(0:Nw/2)'/Nw;
for k = 1:K
    temp = fft(x(ind:ind+Nw-1).*conj(g));
    c(:,k) = temp(1:Nw/2+1).*exp(-wi*(k-1)*R);
    ind = ind + R;
end

if nargout > 1
    E = ones(Nw/2+1,K);
    for k = 2:K
        E(:,k) = exp(-wi*(k-1)*R);
    end
end