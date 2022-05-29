
function computeDMD_oneModeAnimation(rpm,tsignal,save_directory,name,average_degree,Pulse_blade,T)
    dt=T/length(rpm);
    for kk=1:length(rpm)
        len(kk)=length(rpm{1,kk});
    end
    [q1,q2]=min(len);
    for kk=1:length(rpm)
        VORTALL(:,kk)=reshape(rpm{1,kk}(1:min(len),:),1,[])';%按列
    end

    X = VORTALL(:,1:end-1);
    X1 = VORTALL(:,2:end);
    [U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
    r =3 %length(rpm)-1;  %; truncate at 21 modes
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);
    Atilde = U'*X1*V*inv(S);
    [W,eigs] = eig(Atilde);
    Phi =abs(X1*V*inv(S)*W);
    lambda = diag(eigs); % discrete -time eigenvalues
    omega = log(lambda)/dt; % continuous-time eigenvalues
    abs(omega)
%%  Plot DMD spectrum
    h=figure
    theta = (0:1:100)*2*pi/100;
    plot(cos(theta),sin(theta),'k--') % plot unit circle
    hold on, grid on
    scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
    axis([-1.1 1.1 -1.1 1.1]);
    saveas(h,[save_directory,'\',name,'-d',num2str(average_degree),'r',num2str(r),'.png'])
%% Compute DMD mode amplitudes b
    x = X(:, 1);
    b = pinv(Phi)*x;
    
    
%% DMD reconstruction全模态 及其 单个模态（除去背景）分析-并输出tecplot
     mm1 = size(X, 2); % mm1 = m - 1
     t = (0:mm1 -1)*dt; % time vector
     Xdmd_recon = Phi*diag(b)*exp(omega*t);
    
     %选取前10%
     n=round(length(omega)*0.1);
     if mod(n,2) == 0
         n=n+1;
     end
     Omega=abs(omega); oder_Omeag=sort(Omega(:));
     [m1]=find(Omega<=oder_Omeag(n));
     [m2]=find(Omega>oder_Omeag(n));
     Xdmd_background=Phi(:,m1)*diag(b(m1))*exp(omega(m1)*t);
     Xdmd_left=Phi(:,m2)*diag(b(m2))*exp(omega(m2)*t);
     
     Tsignal.Nvar= 8;
     Tsignal.varnames= tsignal.varnames;
    
        Tsignal.varnames{1, 4}=['Xdmd_recon'];
        Tsignal.varnames{1, 5}=['Xdmd_background'];
        Tsignal.varnames{1, 6}=['Xdmd_left'];
        Tsignal.varnames{1, 7}=['X_real'];
        Tsignal.varnames{1, 8}=['X_real-Xdmd_background'];
        for kk=1:length(t)
            Tsignal.surfaces(kk).zonename= tsignal.surfaces(q2(1)).zonename;
            Tsignal.surfaces(kk).x= tsignal.surfaces(q2(1)).x;
            Tsignal.surfaces(kk).y= tsignal.surfaces(q2(1)).y;
            Tsignal.surfaces(kk).z= tsignal.surfaces(q2(1)).z;  
            Tsignal.surfaces(kk).v(1,:,:)= reshape(abs(Xdmd_recon(:,kk)),min(len),10); 
            Tsignal.surfaces(kk).v(2,:,:)= reshape(abs(Xdmd_background(:,kk)),min(len),10); 
            Tsignal.surfaces(kk).v(3,:,:)= reshape(abs(Xdmd_left(:,kk)),min(len),10); 
            Tsignal.surfaces(kk).v(4,:,:)= reshape(tsignal.surfaces(kk).v(1,:,:),min(len),10); 
            Tsignal.surfaces(kk).v(5,:,:)= reshape(tsignal.surfaces(kk).v(1,:,:),min(len),10)-reshape(abs(Xdmd_background(:,kk)),min(len),10); 
            
            Tsignal.surfaces(kk).solutiontime=kk;
        end        
        mat2tecplot(Tsignal,[save_directory,'\',name,'DMD-refine-ver4-d',num2str(average_degree),'.plt']);
end