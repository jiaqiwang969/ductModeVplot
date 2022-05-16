function r=correlation_instant(x,y,Nt,max_tau)
 
x=x-mean(x);y=y-mean(y);

Nt=(length(x)-2*max_tau)/2;
r=zeros(2*max_tau+1,Nt);

t=(1:Nt); 
for k=-max_tau:max_tau   
   r(k+max_tau+1,:)=x(t+max_tau).*conj(y(t+k+max_tau));                                  
end

