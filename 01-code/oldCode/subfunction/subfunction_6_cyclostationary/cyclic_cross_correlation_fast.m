function R=cyclic_cross_correlation_fast(x,y,T,max_tau)

if nargin~=4
  error('Incorrect number of arguments for function cyclic_cross_correlation_fast');
end
if T(1)<1
  error('Synchronous period must be larger than 1 in function cyclic_cross_correlation_fast');
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

if length(T)==1
  if T~=floor(T)
    cp=1;np=1;
    while cp+T<length(x)
      cp=cp+floor(T);
      np=np+T;
      if (np-cp)>1
	x=[x(1:cp-1) x(cp+1:length(x))];
	y=[y(1:cp-1) y(cp+1:length(y))];
	np=np-1;
      end
    end
  end
  n=floor((length(x)-2*max_tau-1)/T);
else
  rot_positions=T;
  T=floor(median(diff(rot_positions)));
  nx=[];
  ny=[];
  n=length(rot_positions)-2;
  for k=1:n;
    cp=rot_positions(k);
    nx=[nx  x(cp:cp+T-1)];
    ny=[ny  y(cp:cp+T-1)];
  end
  nx=[nx x(rot_positions(n+1):rot_positions(n+1)+2*max_tau+1)];
  x=nx;
  ny=[ny y(rot_positions(n+1):rot_positions(n+1)+2*max_tau+1)];
  y=ny;
end

x=x-mean(x);y=y-mean(y);             


r=zeros(2*max_tau+1,floor(T));
temp=zeros(floor(T),n);
t=(1:floor(T)*n); 
for k=-max_tau:max_tau   
   temp(:)=x(t+max_tau).*y(t+k+max_tau);
   [rowt,colt]=size(temp);
   
   if colt==1                               
       r(k+1+max_tau,:)=temp';
   else
       r(k+1+max_tau,:)=mean(temp');        
  end                                         
                                          
   
end


 R=zeros(2*max_tau+1,floor(T));

for k=-max_tau:max_tau
    R(k+1+max_tau,:)=exp(-i*pi*((0:floor(T)-1)/T)*k).*fft(r(k+1+max_tau,:)-mean(r(k+1+max_tau,:)))/T;   
end

R=R(:,1:floor(floor(T)/2));

if nargout==0  
    

  [r,c]=size(R);
  dispR=R(1:r,1:floor(c));
  dispR=fftshift(dispR);
  xl=(1:floor(c))-floor(c/2)-1;
  yl=(1:r)-1;

  xl=xl/c;yl=yl/r;
  figure(gcf+1);
  contour(xl,yl,abs(dispR));
  title('Spectral Correlation Density');
  xlabel('alpha * pi radians');
  ylabel('tao ');
end



