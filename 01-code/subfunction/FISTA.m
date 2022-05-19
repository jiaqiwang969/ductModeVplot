function [R_matrix, para_val]  = FISTA(Li,N_iter,D_measured,SC,mu_final,mu,psi_B);

% Fast iterative shrinkage thresholding algorithm (FISTA)
% (Yu liang prepared in 09/02/2015) 

% input parameters:
% Li is the step size Li = 1.2 or 0.5
% N_iter is maximum iteration for each regularization step N_iter = 10 or 20
% D_measured is measured spectral matrix
% SC stopping criterion 
% mu is initial regularization parameter mu = norm(D_measured,'fro')
% mu_final is final regularization parameter mu_final = 1e-16*mu
% psi_B is the constructed Fourier basis  

D_datamissing = D_measured; % take the measurements data which is a data missing matrix
D_measured(find(D_measured~=0)) = 1;
PRDmatrix = D_measured;  % the positions which the measurements are nonzeros

% initial parameters 
Xk = zeros(size(D_measured,1),size(D_measured,2)); % completed matrix
Yk = Xk; % 
tk = 1;
outer_iter = 1;
psiB_T = psi_B';

% start the iterations
tic
while mu > mu_final
    mu = max(mu * 0.7, mu_final);
for i = 1:N_iter
    Gk = Yk + Li.* (D_datamissing - Yk.*PRDmatrix); % gradient descent 
    
    Gk = psi_B*Gk*psiB_T; % smooth the spectral matrix
    
    % low rank estimation by SVD  
    % Gk = (Gk + Gk')/2;
    % [U,S,V] = svd(Gk);
    % Xkk = U*diag(max(diag(S) - mu*Li,0))*V'; 
    
    % low rank estimation by EVD 
    Gk = (Gk + Gk')/2;
    [V,S]= eig(Gk);
    [Y,I]= sort(diag(S),'descend'); 
    S = diag(Y);
    V = V(:,I);
    Xkk = V*diag(max(diag(S) - mu*Li,0))*V';

    tk1 = 0.5 + 0.5*sqrt(1+4*tk^2);
    Yk = Xkk + ((tk - 1)/tk1)*(Xkk - Xk); % low rank matrix update
    Xk = Xkk;
    tk = tk1; 
end

    err = norm(D_datamissing - Yk.*PRDmatrix,'fro')/norm(D_datamissing,'fro');
    if  err < SC;
        break
    end
    
    total_iter = N_iter * outer_iter; 
    outer_iter = outer_iter + 1;
end

R_matrix = Xk; % complted spectral matrix
para_val.cost_time = toc; % cost time
para_val.mu = mu; % final mu
para_val.err = err; % residual error
para_val.total_iter = total_iter; % total iteration steps
end