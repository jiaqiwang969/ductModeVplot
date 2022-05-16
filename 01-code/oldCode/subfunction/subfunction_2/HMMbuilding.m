 function [LL_r, prior1_r, transmat1_r, mu1_r, Sigma1_r, mixmat1_r]=HMMbuilding(Q,M,traindata,IterationNo)
traindata=cell2mat(traindata);
 [O,T,nex]=size(traindata);
%********************************HMMѵ��*******************************************************************
% Q = 3;   %״̬��ĿN
% M = 2;   %�����M   N*M<T*nex

cov_type = 'full';
%**************ѵ���������***********************************%

%*******�����ʼ��   ��ʼ�ֲ�pi �� ״̬ת�Ƹ��ʾ���A************%
prior0_r = -(rand(Q,1));    %ά��ΪQ��1���� N��1
transmat0_r = mk_stochastic(rand(Q,Q));%ά��ΪQ��Q���� N��N %ÿ�м�����Ϊ1
%ѵ������Ϊ  train_data{1}
%***����k-mean�㷨��ʼ��CHMMģ�Ͳ��� ��ֵmu��Э����U��sigma�����������C��mixmat��*******%
[mu0_r, Sigma0_r] = mixgauss_init(Q*M, traindata, cov_type);
mu0_r = reshape(mu0_r, [O Q M]);
Sigma0_r = reshape(Sigma0_r, [O O Q M]);
mixmat0_r = mk_stochastic(rand(Q,M));
% improve guess of parameters using EM
[LL_r, prior1_r, transmat1_r, mu1_r, Sigma1_r, mixmat1_r] = ...
    mhmm_em(traindata, prior0_r, transmat0_r, mu0_r, Sigma0_r, mixmat0_r, 'max_iter', IterationNo);
%*************************************************************************%

