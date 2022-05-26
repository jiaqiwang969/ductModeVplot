function [x,gg] = irSTFT(c,g,R,I)
% [x,gg] = irSTFT(c,g,R,I)
% Gabor���任
% ���������
% c �� Gabor�任������ϵ������
% g �� ʱ��Ĵ��ں���
% R �Ǵ��ں����ƶ��ľ��루��Ӧ�����任��
% ��������������3���ź�xֻ�ؽ���Ӧ�� I ��Ƶ�ʵ���ź�
% ���������gg ���źű߽紦�ضϵ���Ŀ

[Nw,K] = size(c);
Nw = 2*(Nw-1);
g = g(:);
g = g/sum(g(1:R:Nw).^2);
L = (K-1)*R + Nw;

if nargin < 4
    c = ifft([c;conj(c(Nw/2:-1:2,:))]);
else
    if isempty(I)
        c = ifft([c;conj(c(Nw/2:-1:2,:))]);
    else
        C = zeros(Nw,K);
        n = [0:Nw-1]';
        for k = I
            C = C + exp(2i*pi*(k-1)*n/Nw)*c(k,:)/Nw;
        end
        c = C;
        clear C
    end
end

x = zeros(L,1);
ind = 1:Nw;
for k = 1:K
    x(ind) = x(ind) + g.*c(mod(ind-1,Nw)+1,k);
    ind = ind + R;
end
x = real(x);

if nargout > 1
    gg = zeros(L,1);
    ind = 1:Nw;
    for k = 1:K
        gg(ind) = gg(ind) + g.^2;
        ind = ind + R;
    end
end
