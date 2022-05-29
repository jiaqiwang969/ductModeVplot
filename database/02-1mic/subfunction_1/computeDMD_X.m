
function computeDMD_X(fs,signal_X,modeNumber)
[U, S, V] = svd(signal_X(1:50000,1));
pc1 = U(:, 1); % first PCA mode
pc2 = U(:, 2); % second PCA mode
time_pc1 = V(:, 1); % temporal evolution of pc1
time_pc2 = V(:, 2); % temporal evolution of pc2

figure 
plot(pc1)
figure 
plot(signal_X(1:50000,1))
dt=1/fs;
VORTALL=signal_X';
X = VORTALL(:,1:end-1);
X1 = VORTALL(:,2:end);

[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)

U = U(:,1:modeNumber);
S = S(1:modeNumber,1:modeNumber);
V = V(:,1:modeNumber);
Atilde = U'*X1*V*inv(S);
[W,eigs] = eig(Atilde);
Phi =real(X1*V*inv(S)*W);

lambda = diag(eigs); % discrete -time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues

%%  Plot DMD spectrum
h=figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

%% Compute DMD mode amplitudes b
x = X(:, 1);
b = Phi\x;
%% DMD reconstruction
mm1 = size(X, 2); % mm1 = m - 1
time_dynamics = zeros(modeNumber, mm1);
t = (0:mm1 -1)*dt; % time vector
for iter = 1:length(t)
time_dynamics (:,iter )=(b.*exp(omega*t(iter )));
end
Xdmd = Phi * time_dynamics ;

%单个模态分析
k=1;
for iter = 1:length(t)
Xdmd_k(:,iter )=Phi(:,k)*exp(omega(k)*t(iter));
end



end