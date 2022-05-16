function x=createsignal(fs,fa,thitaA,ampA,fb,thitaB,ampB,ampC,N)
if nargin< 8 
    error('There must have eight input elements.');
end

[row,col]=size(fb);
if row>col
    fb=fb';
end
[row,col]=size(thitaB);
if row>col
    thitaB=thitaB';
end
[row,col]=size(ampB);
if row>col
    ampB=ampB';
end

[rowfb,colfb]=size(fb);
[rowB,colB]=size(thitaB);
[rowampB,colampB]=size(ampB);

if colfb>colB
    dnum=colfb-colB;
    for i=1:dnum
        thitaB=[thitaB,0];
    end
end
if colfb>colampB
    dnum=colfb-colampB;
    for i=1:dnum
        ampB=[ampB,1];
    end
end    

dt=1/fs;
n=0:N-1;
t=n*dt;

sB=zeros(colfb,N);
sumB=zeros(1,N);
for i=1:colfb
    sB(i,:)=ampB(i)*cos(2*pi*fb(i)*t+thitaB(i));
end
if colfb==1
    sumB=sB(1,:);
else
    sumB=sum(sB);
end


width_f=max(fb);
band_pass=[fa-width_f fa+width_f]/fs;
x=ampA*(ampC+sumB).*cos(2*pi*fa*t+thitaA);

if nargout==0
    figure(gcf+1);
    plot(t,x);
    xlabel('time[s]');
    ylabel('amplitude');
    title('signal');
end