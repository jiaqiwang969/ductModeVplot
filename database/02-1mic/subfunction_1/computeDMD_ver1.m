
function computeDMD_ver1(rpm,tsignal,save_directory,name,average_degree,Pulse_blade,T4oneRound)
dt=T4oneRound/length(rpm);
for kk=1:length(rpm)
    len(kk)=length(rpm{1,kk});
end
[q1,q2]=min(len)
for kk=1:length(rpm)
    VORTALL(:,kk)=reshape(rpm{1,kk}(1:min(len),:),1,[])';%按列
end

X = VORTALL(:,1:end-1);
X2 = VORTALL(:,2:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r =5;%length(rpm)-1;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi =abs(X2*V*inv(S)*W);

lambda = diag(eigs); % discrete -time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
%% Compute DMD mode amplitudes b
x = X(:, 1);
b = Phi\x;
%% DMD reconstruction
mm1 = size(X, 2); % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:10*mm1 -1)*dt/10; % time vector
for iter = 1:length(t)
time_dynamics (:,iter )=(b.*exp(omega*t(iter )));
end
Xdmd = Phi * time_dynamics ;

%单个模态分析
k=4;
for iter = 1:length(t)
Xdmd_k(:,iter )=Phi(:,k)*exp(omega(k)*t(iter));
end





tsignal2.Nvar= tsignal.Nvar;
tsignal2.varnames= tsignal.varnames;
% for kk=1:r
for kk=1:iter
tsignal2.surfaces(kk).zonename= tsignal.surfaces(q2(1)).zonename;
tsignal2.surfaces(kk).x= tsignal.surfaces(q2(1)).x;
tsignal2.surfaces(kk).y= tsignal.surfaces(q2(1)).y;
tsignal2.surfaces(kk).z= tsignal.surfaces(q2(1)).z;
% tsignal2.surfaces(kk).v= reshape(Xdmd(:,kk),1,min(len),10); 
% tsignal2.surfaces(kk).v(1,Pulse_blade,:)=1.1*max(Phi(:,kk));
tsignal2.surfaces(kk).v= reshape(abs(Xdmd_k(:,kk)),1,min(len),10); 
%tsignal2.surfaces(kk).v(1,Pulse_blade,:)=1.1*max(real(Xdmd(:,kk)));
tsignal2.surfaces(kk).solutiontime=kk;
end
mat2tecplot(tsignal2,[save_directory,'\',name,'-d',num2str(average_degree),'r',num2str(r),'-reconstruction','.plt']);

%%  Plot DMD spectrum
h=figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);
saveas(h,[save_directory,'\',name,'-d',num2str(average_degree),'r',num2str(r),'.png'])
end