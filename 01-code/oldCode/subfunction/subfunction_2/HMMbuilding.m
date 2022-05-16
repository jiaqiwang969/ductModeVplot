 function [LL_r, prior1_r, transmat1_r, mu1_r, Sigma1_r, mixmat1_r]=HMMbuilding(Q,M,traindata,IterationNo)
traindata=cell2mat(traindata);
 [O,T,nex]=size(traindata);
%********************************HMM训练*******************************************************************
% Q = 3;   %状态数目N
% M = 2;   %混合数M   N*M<T*nex

cov_type = 'full';
%**************训练正常轴承***********************************%

%*******随机初始化   初始分布pi 和 状态转移概率矩阵A************%
prior0_r = -(rand(Q,1));    %维数为Q×1，即 N×1
transmat0_r = mk_stochastic(rand(Q,Q));%维数为Q×Q，即 N×N %每行加起来为1
%训练数据为  train_data{1}
%***采用k-mean算法初始化CHMM模型参数 均值mu、协方差U（sigma）、增益矩阵C（mixmat）*******%
[mu0_r, Sigma0_r] = mixgauss_init(Q*M, traindata, cov_type);
mu0_r = reshape(mu0_r, [O Q M]);
Sigma0_r = reshape(Sigma0_r, [O O Q M]);
mixmat0_r = mk_stochastic(rand(Q,M));
% improve guess of parameters using EM
[LL_r, prior1_r, transmat1_r, mu1_r, Sigma1_r, mixmat1_r] = ...
    mhmm_em(traindata, prior0_r, transmat0_r, mu0_r, Sigma0_r, mixmat0_r, 'max_iter', IterationNo);
%*************************************************************************%

