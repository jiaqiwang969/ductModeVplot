% sparse_in_Spatial domain.m  by wjq 2018-06-13
  
%  
%This code demonstrate compressive sensing example. In this  
%example the signal is sparse in Circumferential mode and random samples  
%are taken in Spatial domain.  
function [rotor_speed]=mode_dB_universal_CS(K,signal,fs,Scale_data,object,nk,objectName,testTime,fname)
    the_freq = [0:length(signal)/2.56 - 1]*fs/length(signal);  %数据频域离散刻度
    temp_freq =fft(signal(:,:))*2./length(signal);
    temp_freq=temp_freq(1:length(signal)/2.56,1:nk);
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
    Fk=[1/29,1/2,1,2,3];
    for k=1:5
        xf=max((temp_freq(floor(rotor_speed/60*29*Fk(k)/(fs/length(signal)))+[floor(-2/(fs/length(signal))):floor(2/(fs/length(signal)))],:)))';

        %creating dft matrix   
        m = [-15:16].';                 % Angular increment.
        theta_k = (0:2*pi/nk:2*pi-2*pi/nk).* 1i;  % Column.
        B = exp(-m*theta_k)./nk;              % Exponentiation of outer product.
        Binv=inv(B);  % The inverse discrete Fourier transform matrix, Binv, equals CONJ(dftmtx(N))/N.   
        %Selecting random rows of the DFT matrix  
        q=randperm(nk);  
        %creating measurement matrix  
        A=Binv(q(1:K),:);    % 在IDFT矩阵中任选K=256行  
        %taking random time measurements  
        %y=(A*xf);   % 对x的fft后的xf(1024-by-1)的数据做IDFT得到256个时域稀疏采样值，通过plot(real(y))和原来的x对比，注意如何在时域中取K=256个采样值  
        y=xf(q(1:K),1);
        %Calculating Initial guess  
        x0=A.'*y;  % 注意：待恢复时域信号xprec的DFT值xp的估计初值x0如何给？ y 是时域稀疏采样值  
        %Running the recovery Algorithm 
        tic  
        xp(:,k)=l1eq_pd(x0,A,[],y,1e-5); %恢复的xp是频域信号  
        toc  
    end
bar(m,abs(xp));hold on
legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',16);
set(gca,'XTick',m);
set(gca,'Ygrid','on') 
title({[testTime,'-模态分析'];[num2str(K),'支传感器-',fname,'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
   
end



%%%%%%%%%%%%%%%%%%飘逸的分割线%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
% 用到的函数  
% l1eq_pd.m  
%  
% Solve  
% min_x ||x||_1  s.t.  Ax = b  
%  
% Recast as linear program  
% min_{x,u} sum(u)  s.t.  -u <= x <= u,  Ax=b  
% and use primal-dual interior point method  
%  
% Usage: xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)  
%  
% x0 - Nx1 vector, initial point.  
%  
% A - Either a handle to a function that takes a N vector and returns a K   
%     vector , or a KxN matrix.  If A is a function handle, the algorithm  
%     operates in "largescale" mode, solving the Newton systems via the  
%     Conjugate Gradients algorithm.  
%  
% At - Handle to a function that takes a K vector and returns an N vector.  
%      If A is a KxN matrix, At is ignored.  
%  
% b - Kx1 vector of observations.  
%  
% pdtol - Tolerance for primal-dual algorithm (algorithm terminates if  
%     the duality gap is less than pdtol).    
%     Default = 1e-3.  
%  
% pdmaxiter - Maximum number of primal-dual iterations.    
%     Default = 50.  
%  
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.  
%     Default = 1e-8.  
%  
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored  
%     if A is a matrix.  
%     Default = 200.  
%  
% Written by: Justin Romberg, Caltech  
% Email: jrom@acm.caltech.edu  
% Created: October 2005  
%  
  
function xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)  
  
largescale = isa(A,'function_handle');  
  
if (nargin < 5), pdtol = 1e-3;  end  
if (nargin < 6), pdmaxiter = 50;  end  
if (nargin < 7), cgtol = 1e-8;  end  
if (nargin < 8), cgmaxiter = 200;  end  
  
N = length(x0);  
  
alpha = 0.01;  
beta = 0.5;  
mu = 10;  
  
gradf0 = [zeros(N,1); ones(N,1)];  
  
% starting point --- make sure that it is feasible  
if (largescale)  
  if (norm(A(x0)-b)/norm(b) > cgtol)  
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');  
    AAt = @(z) A(At(z));  
    [w, cgres, cgiter] = cgsolve(AAt, b, cgtol, cgmaxiter, 0);  
    if (cgres > 1/2)  
      disp('A*At is ill-conditioned: cannot find starting point');  
      xp = x0;  
      return;  
    end  
    x0 = At(w);  
  end  
else  
  if (norm(A*x0-b)/norm(b) > cgtol)  
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');  
    opts.POSDEF = true; opts.SYM = true;  
    [w, hcond] = linsolve(A*A', b, opts);  
    if (hcond < 1e-14)  
      disp('A*At is ill-conditioned: cannot find starting point');  
      xp = x0;  
      return;  
    end  
    x0 = A'*w;  
  end    
end  
x = x0;  
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));  
  
% set up for the first iteration  
fu1 = x - u;  
fu2 = -x - u;  
lamu1 = -1./fu1;  
lamu2 = -1./fu2;  
if (largescale)  
  v = -A(lamu1-lamu2);  
  Atv = At(v);  
  rpri = A(x) - b;  
else  
  v = -A*(lamu1-lamu2);  
  Atv = A'*v;  
  rpri = A*x - b;  
end  
  
sdg = -(fu1'*lamu1 + fu2'*lamu2);  
tau = mu*2*N/sdg;  
  
rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);  
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];  
resnorm = norm([rdual; rcent; rpri]);  
  
pditer = 0;  
done = (sdg < pdtol) | (pditer >= pdmaxiter);  
while (~done)  
    
  pditer = pditer + 1;  
    
  w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv;  
  w2 = -1 - 1/tau*(1./fu1 + 1./fu2);  
  w3 = -rpri;  
    
  sig1 = -lamu1./fu1 - lamu2./fu2;  
  sig2 = lamu1./fu1 - lamu2./fu2;  
  sigx = sig1 - sig2.^2./sig1;  
    
  if (largescale)  
    w1p = w3 - A(w1./sigx - w2.*sig2./(sigx.*sig1));  
    h11pfun = @(z) -A(1./sigx.*At(z));  
    [dv, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);  
    if (cgres > 1/2)  
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');  
      xp = x;  
      return  
    end  
    dx = (w1 - w2.*sig2./sig1 - At(dv))./sigx;  
    Adx = A(dx);  
    Atdv = At(dv);  
  else  
    w1p = -(w3 - A*(w1./sigx - w2.*sig2./(sigx.*sig1)));  
    H11p = A*(sparse(diag(1./sigx))*A');  
    opts.POSDEF = true; opts.SYM = true;  
    [dv,hcond] = linsolve(H11p, w1p);  
    if (hcond < 1e-14)  
      disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');  
      xp = x;  
      return  
    end  
    dx = (w1 - w2.*sig2./sig1 - A'*dv)./sigx;  
    Adx = A*dx;  
    Atdv = A'*dv;  
  end  
    
  du = (w2 - sig2.*dx)./sig1;  
    
  dlamu1 = (lamu1./fu1).*(-dx+du) - lamu1 - (1/tau)*1./fu1;  
  dlamu2 = (lamu2./fu2).*(dx+du) - lamu2 - 1/tau*1./fu2;  
    
  % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0  
  indp = find(dlamu1 < 0);  indn = find(dlamu2 < 0);  
  s = min([1; -lamu1(indp)./dlamu1(indp); -lamu2(indn)./dlamu2(indn)]);  
  indp = find((dx-du) > 0);  indn = find((-dx-du) > 0);  
  s = (0.99)*min([s; -fu1(indp)./(dx(indp)-du(indp)); -fu2(indn)./(-dx(indn)-du(indn))]);  
    
  % backtracking line search  
  suffdec = 0;  
  backiter = 0;  
  while (~suffdec)  
    xp = x + s*dx;  up = u + s*du;   
    vp = v + s*dv;  Atvp = Atv + s*Atdv;   
    lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;  
    fu1p = xp - up;  fu2p = -xp - up;    
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];  
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);  
    rpp = rpri + s*Adx;  
    suffdec = (norm([rdp; rcp; rpp]) <= (1-alpha*s)*resnorm);  
    s = beta*s;  
    backiter = backiter + 1;  
    if (backiter > 32)  
      disp('Stuck backtracking, returning last iterate.  (See Section 4 of notes for more information.)')  
      xp = x;  
      return  
    end  
  end  
    
    
  % next iteration  
  x = xp;  u = up;  
  v = vp;  Atv = Atvp;   
  lamu1 = lamu1p;  lamu2 = lamu2p;  
  fu1 = fu1p;  fu2 = fu2p;  
    
  % surrogate duality gap  
  sdg = -(fu1'*lamu1 + fu2'*lamu2);  
  tau = mu*2*N/sdg;  
  rpri = rpp;  
  rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);  
  rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];  
  resnorm = norm([rdual; rcent; rpri]);  
    
  done = (sdg < pdtol) | (pditer >= pdmaxiter);  
    
  disp(sprintf('Iteration = %d, tau = %8.3e, Primal = %8.3e, PDGap = %8.3e, Dual res = %8.3e, Primal res = %8.3e',...  
    pditer, tau, sum(u), sdg, norm(rdual), norm(rpri)));  
  if (largescale)  
    disp(sprintf(' CG Res = %8.3e, CG Iter = %d', cgres, cgiter));  
  else  
    disp(sprintf(' H11p condition number = %8.3e', hcond));  
  end  
    
end  
end