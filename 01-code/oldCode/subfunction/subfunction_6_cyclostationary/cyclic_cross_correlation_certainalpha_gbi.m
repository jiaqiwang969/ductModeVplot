function R=cyclic_cross_correlation_certainalpha_gbi(x,y,alpha,fs,Nt,max_tau)

if nargin~=6
  error('Incorrect number of arguments for function cyclic_cross_correlation_fast');
end

if length(x)>length(y)
  x=x(1:length(y));
end
[rows,cols]=size(x);
if rows>cols
  x=x';
end
[rows,cols]=size(y);
if rows>cols
  y=y';
end
 nalpha=floor(alpha*Nt/fs);

n=floor((length(x)-2*max_tau)/Nt);

x=x-mean(x);y=y-mean(y);

r=zeros(2*max_tau+1,Nt);
temp=zeros(Nt,n);
t=(1:Nt*n); 
for k=-max_tau:max_tau   
   temp(:)=x(t+max_tau).*conj(y(t+k+max_tau));
   
   if n==1                  
       r(k+1+max_tau,:)=temp';
   else
       r(k+1+max_tau,:)=mean(temp');       
                                                
   end                                         
end

 R=zeros(2*max_tau+1,1);

taul=-max_tau:max_tau;
factortmp=exp(-i*pi*nalpha/Nt*taul);
 R=factortmp.*(exp(-i*pi*2*nalpha*(0:Nt-1)/Nt)*r'/Nt); 

if nargout==0  
    drawresult(abs(R),'R','contour',1);
end



