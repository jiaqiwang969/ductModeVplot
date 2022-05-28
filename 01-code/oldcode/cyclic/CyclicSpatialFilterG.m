function [G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,epsilon)
% [G,Da] = CyclicSpatialFilterG(Syya,Syaya,Ns,epsilon)
% 构建循环维纳滤波器G, 从而提取循环平稳信号
% 输入:
% Syya = 循环互功率谱矩阵, (y 和 ya 之间的互谱矩阵)
% Syaya = 循环互功率谱矩阵, (ya 和 ya 之间的互谱矩阵)
% Ns = 循环平稳声源的个数
% epsilon (optional) = 矩阵求逆时, 做伪逆的动态范围, 缺省为1e-12

% 输出:
% G = 分离滤波器
% Da = 矩阵 Syya*inv(Syaya)*Syya' 的(对应不同频率的)特征值向量

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