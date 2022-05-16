function drawresult(RS,ShowRS,Drtype,fs,linespace)

if nargin==1
    ShowRS='R'; 
    Drtype='contour';
end 
if nargin<5
    linespace=1;
end

[r,c]=size(RS);     
  dispRS=RS(1:r,1:c);
  
  xl=(1:c+1)-1;              
  xl=xl*fs/(2*c);         
 yl=(1:r+1)-1;  
  if ShowRS=='R'
      yl=yl-floor(r/2); 
      yl=yl/fs;
      titletxt='cyclic Correlation';
      xlabeltxt='cyclic frequency[Hz]';
      ylabeltxt='time delays[s]';
  elseif ShowRS=='S'
         yl=yl*fs/(2*r);
      titletxt='Spectral Correlation Density';
      xlabeltxt='cyclic frequency[Hz]';
      ylabeltxt='frequency[Hz]';
      zlabeltxt='CSD';

  else
  end
   
  
         dispRS=[dispRS' zeros(c,1)];
       dispRS=[dispRS' zeros(r+1,1)]; 
       
   
        tmp=dispRS;
       tmp(abs(tmp)<0.00005)=NaN;
  figure(gcf+1);
  if strcmp(Drtype,'contour')
contour(xl,yl,abs(dispRS));
axis tight;
    set(gcf,'color','white');
  elseif strcmp(Drtype,'surf')
      surf(xl,yl,abs(dispRS));
      shading flat;
 axis tight;
    set(gcf,'color','white');

  elseif strcmp(Drtype,'waterfall')
      waterfall(xl,yl,abs(dispRS));
 elseif strcmp(Drtype,'plot3')
     [X,Y]=meshgrid(0:c-1+1,0:r-1+1);           
     plot3(X,Y,abs(tmp));
     axis tight;
    set(gcf,'color','white');
 elseif strcmp(Drtype,'mesh')
      mesh(xl,yl,abs(dispRS));
      shading flat;
      axis tight;
    set(gcf,'color','white');
  else
  end
 title(titletxt);
 xlabel(xlabeltxt);
 ylabel(ylabeltxt);
zlabel('amplitude');