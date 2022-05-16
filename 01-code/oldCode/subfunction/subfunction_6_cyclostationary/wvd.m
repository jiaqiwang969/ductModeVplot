function W=wvd(S)

W=fft(S')';
[r,c]=size(W);

if nargout==0  
  dispW=(W(1:(r+1)/2,:));
  x=(1:c)-1;
  y=(1:(r+1)/2)-1;
  y=y/r;
  contour(x,y,abs(dispW))
  title('Wigner-Ville Time-Frequency Distribution')
  xlabel('Time/Samples')
  ylabel('Frequency * pi radians')
end


