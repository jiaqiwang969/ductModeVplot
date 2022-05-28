function [x,cx] = CyclicSpatialFiltering(y,G,phi,alpha,orders,Window,R)
% [x,cx] = CyclicSpatialFiltering(y,G,phi,alpha,orders,Window,R)
% ���ڶ�ʱ����Ҷ�任�ķ�ʽ��ѭ��ƽ���źŽ�����ȡ
% 
% ����:
% y =  �����ź�, ����ÿһ��Ϊһ��ͨ���ź�;
% G = �����˲���
% phi = ˲ʱ��λ����λ�ǻ���
% alpha = ˲ʱƵ�� (��λΪHz)����ͨ�����Բ���Ƶ�ʹ�һ�� = (d phi(t)/dt)/(2 pi*Fs)
% orders = ѭ��Ƶ�ʵĽ���
% ע��: ������������ phi �� alpha �Ǿ���, ˵��ÿһ�ж�Ӧ��ͬ��ѭ����λ��ѭ��Ƶ��;
% ����, alpha �Ľ���ֱ�Ӵ��������� 'orders'��ѡȡ
% Window (Gabor �任����) 
% = Welch ������ʹ�õĴ��ڳ���, ��������Gabor�任��ϵ��,ȱʡΪHanning(Nw)
% R (��ѡ,Gabor �任����) = �����ƶ�����, default = Nw/4, ��Ӧ�ظ��Ĵ���Ϊ3/4*Nw
% �ú������� 'CSM_W', 'rgt'  'irgt' �Ӻ���.
%
% ���:
% --------
% x = �˲�����ź�, ����ȡ��ѭ��ƽ�ȳɷ�
% cx = ѭ��ƽ�ȳɷ�Gabor��ϵ��

[L,M] = size(y);
N = size(G,1);

if nargin < 6
    Nw = size(G,3);
    Nw = 2*(Nw-1);
    R = fix(Nw/4);
    Window = hanning(Nw);
else
    Nw = length(Window);
end

if min(size(phi)) > 1
    Harmonics = 0;
    Na = size(phi,2);
else
    Harmonics = 1;
    Na = length(orders);
end

K = floor((L-Nw)/R+1);
cx = zeros(N,Nw/2+1,K);
for ka = 1:Na
    if Harmonics == 1
        temp = exp(orders(ka)*1i*phi(:)).*alpha(:);
    else
        temp = exp(1i*phi(:,ka)).*alpha(:,ka);
    end
    for m = 1:M
        cya = rgt(y(:,m).*temp,Window,R);
        for n = 1:N
            cx(n,:,:) = squeeze(cx(n,:,:)) + repmat(squeeze(G(n,(ka-1)*M+m,:)),1,K).*cya;
        end
    end
end

L = (K-1)*R + Nw;
x = zeros(L,N);
for n = 1:N
    x(:,n) = irgt(squeeze(cx(n,:,:)),Window,R);
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cya = zeros(Na*M,Nw/2+1,K);
for ka = 1:Na
    temp = exp(orders(ka)*1i*phi(:)).*alpha(:);
    for m = 1:M
        cya((ka-1)*M+m,:,:) = rgt(y(:,m).*temp,Window,R);
    end
end

cx = zeros(M,Nw/2+1,K);
for k = 1:Nw/2+1
    cx(:,k,:) = squeeze(G(:,:,k))*squeeze(cya(:,k,:));
end

L = (K-1)*R + Nw;
x = zeros(L,M);
for m = 1:M
    x(:,m) = irgt(squeeze(cx(m,:,:)),Window,R);
end


