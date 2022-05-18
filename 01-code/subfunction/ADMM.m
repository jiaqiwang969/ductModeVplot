function [R_matrix,err]  = ADMM(D_measured,psi_B, SC, mIter, gama, mu, alpha,SNR,Nl)

% ADMM method  
% input parameters:
% mIter is maximum iteration 迭代次数
% D_measured is measured spectral matrix 测量光谱矩阵
% SC stopping criterion  截止标准
% alpha is initial regularization parameter 初始正则化参数
% psi_B is the constructed Fourier basis
% mu is penalty parameter 惩罚参数
% gama is relaxation parameter
% eta in [0,1)

eta=0.05;
D_datamissing = D_measured;
D_measured(find(D_measured~=0)) = 1;
Omega = D_measured; 
Inv_Omega = ones(size(Omega)) - Omega;
[m, n] = size(Omega);
% initial parameters 
Yk{1} =rand(size(D_measured,1),size(D_measured,2));
Xk{1}=Yk{1};
Pi_k{1} =zeros(size(D_measured,1),size(D_measured,2));

psiB_T = psi_B';
it = 1;
NormDMissing = norm(D_datamissing,'fro');

%lambda_max  =   5; 

while it <  mIter+1
    % update Xk, low rank estimation by EVD
     Gk{it} = Yk{it} - (1/mu)*Pi_k{it};
     Gk{it} = psi_B *  Gk{it} * psiB_T; % smooth the spectral matrix 
     Gk{it} = (Gk{it} +Gk{it}')/2; %
     t= alpha/mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Uk,Wk,Vk] = eig(Gk{it});  
    Xk{it+1} = Uk*max(Wk - t, 0)*Vk'; %first step
    Ek = Xk{it+1} + (1/mu)*Pi_k{it};
    Yk{it+1} = 1/(mu+1)*(D_datamissing + mu*Ek).*Omega;  %second step
    Yk{it+1} = Yk{it+1} + Ek.*Inv_Omega;
    Pi_k{it+1} = Pi_k{it} + gama* mu*(Xk{it+1} - Yk{it+1});   %third step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    err = norm(D_datamissing - Xk{it}.*Omega,'fro')/NormDMissing;
    if  err < SC
        break
    end   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  it =it+1;   
end
R_matrix = Xk{it};
end

