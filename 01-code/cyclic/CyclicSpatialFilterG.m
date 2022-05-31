function [G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,epsilon)
% [G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,epsilon)
% ����ѭ��ά���˲���G, �Ӷ���ȡѭ��ƽ���ź�
% ����:
% Syya = ѭ���������׾���, (y �� ya ֮��Ļ��׾���)
% Syaya = ѭ���������׾���, (ya �� ya ֮��Ļ��׾���)
% Ns = ѭ��ƽ����Դ�ĸ���
% epsilon (optional) = ��������ʱ, ��α��Ķ�̬��Χ, ȱʡΪ1e-12

% ���:
% G = �����˲���
% Da = ���� Syya*inv(Syaya)*Syya' ��(��Ӧ��ͬƵ�ʵ�)����ֵ����

if nargin < 4
    epsilon = 1e-12;
end

[M,Ma,Nw2] = size(Syya.S);

G = zeros(M,Ma,Nw2);
Da = zeros(M,Nw2);
for k = 1:Nw2
    [U,D] = eig(squeeze(Syaya.S(:,:,k)));
    [D,I] = sort(diag(D),'descend');
    U = U(:,I);
    I = find(D >= epsilon*D(1));
    invSyaya = U(:,I)*diag(1./D(I))*U(:,I)';
    
    C = squeeze(Syya.S(:,:,k))*invSyaya*squeeze(Syya.S(:,:,k))';
    C = .5*(C + C');
    [U,D] = eig(C);
    [Da(:,k),I] = sort(diag(D),'descend');
    U = U(:,I(1:Ns));
    G(:,:,k) = U*U'*squeeze(Syya.S(:,:,k))*invSyaya;
end