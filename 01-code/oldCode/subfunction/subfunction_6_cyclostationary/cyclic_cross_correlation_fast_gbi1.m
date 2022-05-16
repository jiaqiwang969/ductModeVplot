function R=cyclic_cross_correlation_fast_gbi1(x,y,Nt,max_tau)

if nargin~=4
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

n=floor((length(x)-2*max_tau)/Nt);

x=x-mean(x);y=y-mean(y);

r=zeros(2*max_tau+1,Nt);
temp=zeros(Nt,n);
t=(1:Nt*n); 
for k=-max_tau:max_tau   
   temp(:)=x(t+max_tau).*y(t+k+max_tau);
   [rowt,colt]=size(temp);
   
   if colt==1                               
       r(k+1+max_tau,:)=temp';
   else
       r(k+1+max_tau,:)=mean(temp');     
end                                         
end

 R=zeros(2*max_tau+1,Nt);

for k=-max_tau:max_tau
    R(k+1+max_tau,:)=exp(-i*pi*((0:Nt-1)/Nt)*k).*fft(r(k+1+max_tau,:))/Nt;   
end

if nargout==0  
    drawresult(abs(R),'R','contour');
end



