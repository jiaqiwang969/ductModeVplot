function [Syya,Syaya] = CyclicSpecMat(y,phi,alpha,orders,Window,Noverlap)
% [Syya,Syaya] = CyclicSpecMat(y,phi,alpha,orders,Window,Noverlap)
% ����ѭ���������׾���Syya �� Syaya
% ���� y ��ԭʼ�Ĳ����źţ�ya ��ͨ����ѭ��Ƶ��alpha��չ����չ�ź�

% �������:
% y �������źŵľ�������ÿһ�а���һ��ͨ��������
% phi = ˲ʱ��λ����λ�ǻ���
% alpha = ˲ʱƵ�� (��λΪHz)����ͨ�����Բ���Ƶ�ʹ�һ�� = (d phi(t)/dt)/(2 pi*Fs)
% orders = ѭ��Ƶ�ʵĽ���
% ע��: ������������ phi �� alpha �Ǿ���, ˵��ÿһ�ж�Ӧ��ͬ��ѭ����λ��ѭ��Ƶ��;
% ����, alpha �Ľ���ֱ�Ӵ��������� 'orders'��ѡȡ
% Window = Welch ������ʹ�õĴ��ڳ���
% ������ڳ��ȵ������Ǳ���, �����Ϊ�������ĳ���
% Noverlap (optional) = ��������һ��ѡȡΪ 3/4*Nw, ��3/4�Ĵ��ڳ���

% �������:
% Syya = ѭ���������׾���, (y �� ya ֮��Ļ��׾���)
% Syaya = ѭ���������׾���, (ya �� ya ֮��Ļ��׾���)
% ���� Syya �� Syaya �ǽṹ���͵�����:
%                   Syya.S = ѭ���������׾���
%                   Syya.f = Ƶ������
%                   Syya.K = blocks��ĸ���
%                   Syya.Var_Reduc = �����������
% 
% �ú������� 'CSM_W' �Ӻ���.

[L,M] = size(y);

if length(Window) == 1
    Nw = Window;
    Window = hanning(Nw);
else
    Nw = length(Window);
end
nfft = Nw;
if nargin < 6
    Noverlap = fix(3/4*Nw);
end


if min(size(phi)) > 1
    Harmonics = 0;
    Na = size(phi,2);
else
    Harmonics = 1;
    Na = length(orders);
end

% ya = [];
% for ka = 1:Na
%     ya = [ya y.*repmat(exp(orders(ka)*1i*phi(:)).*alpha(:),1,M)];
% end
ya = zeros(L,Na*M);
for ka = 1:Na
    if Harmonics == 1
         temp = exp(orders(ka)*1i*phi(:)).*alpha(:);
        % temp = exp(orders(ka)*1i*phi(:));
    else
        temp = exp(1i*phi(:,ka)).*alpha(:,ka);
    end
    for m = 1:M
        ya(:,(ka-1)*M+m) = y(:,m).*temp;
    end
end

Syya = CSM_W(y,ya,nfft,Noverlap,Window);
Syaya = CSM_W(ya,ya,nfft,Noverlap,Window);