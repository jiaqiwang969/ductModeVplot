function R=cyclic_cross_correlation(x,y,alpha,max_tau)

if nargin~=4                                                                                    
    
  error('Incorrect number of arguments for function cyclic_cross_correlation');
end
if alpha>2*pi                                                                              
  error('Cyclic frequency must be less than 2 pi in function cyclic_cross_correlation');
end


T=ceil(2*pi/alpha)-1;
lx=length(x);
t=0:lx-1;
R=zeros(max_tau*2+1,T+1);


for tau=-max_tau:2:max_tau
  for k=0:T
    R(tau+1+max_tau,k+1)=mean(x(1:lx-max_tau-tau).*y(max_tau+tau+1:lx) ...
	.*exp(-j*k*alpha*t(1+(max_tau+tau)/2:lx-(max_tau+tau)/2)));
  end
end

t=t+0.5;
for tau=-max_tau+1:2:max_tau
  for k=0:T
    R(tau+1+max_tau,k+1)=mean(x(1:lx-tau-max_tau).*y(max_tau+tau+1:lx) ...
	.*exp(-j*k*alpha*t(1+(max_tau+tau-1)/2:lx-(max_tau+tau+1)/2)));
  end
end



