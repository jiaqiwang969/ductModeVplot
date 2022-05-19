function [R_matrix,para_AP]  = AP_cycle(Hard_Th,D_measured,psi_B,max_it,SC,meanORG,varORG);
                                                             
%                          Cyclic Projection (CP) Algorithm

% input parameter
% Hard_Th: estimated rank of the spectral matrix
% D_measured: measured spectral matrix
% psi_B: reduced spatial basis 
% max_it: maximum iteration steps
% meanORG,varORG: set to be 0; set to be 1 when the random intilization
% case

k = Hard_Th;
psiB_T = psi_B';
PRDmatrix  = find(D_measured~=0); % take the position of the matrix

% method 1: random intilization
% [m,n] = size(D_measured);
% matrixdataX = meanORG + varORG.*(randn(m,n)+1i*randn(m,n));
% matrixdataX = (matrixdataX + matrixdataX')/2;
% matrixdataX(PRDmatrix) = D_measured(PRDmatrix); 

% method 2: smoothing intilization
matrixdataX = psi_B*D_measured*psiB_T ;
matrixdataX(PRDmatrix) = D_measured(PRDmatrix); 

c = 0;
tic
while(c < max_it)
c = c + 1;

% cyclic 1 
matrixdataX = (matrixdataX + matrixdataX')/2;

% low rank estimation by EVD
[V,D]= eig(matrixdataX);
[Y,I]= sort(diag(D),'descend'); 
D = diag(Y);
V = V(:,I);
matrixdataX = V(:,1:k)*D(1:k,1:k)*V(:,1:k)'; % projection to S2

% low rank estimation by SVD
% projection to C2
% [U,D,V] = svd(matrixdataX);
% matrixdataX = U(:,1:k)*D(1:k,1:k)*V(:,1:k)';

% projection to C1
matrixdataX(PRDmatrix) = D_measured(PRDmatrix);

% projection to C3
matrixdataX = psi_B*matrixdataX*psiB_T ; 

err = norm(matrixdataX(PRDmatrix) - D_measured(PRDmatrix),'fro')/norm(D_measured(PRDmatrix),'fro');
if  err < SC % check the stop criteria
    matrixdataX(PRDmatrix) = D_measured(PRDmatrix);
    break
end

end

R_matrix = matrixdataX; % reconstructed matrix
para_AP.tElapsed = toc; % cost time
para_AP.max_it = c; % iteration steps
para_AP.S_valk1 = D(k,k); % kth eigenvalue
para_AP.S_valk2 = D(k+1,k+1); % (k+1)th eigenvalue
para_AP.err = err; % residual error