function [G,Sxdx,Sxx,k]=STFT_LE(x,delay,Nwind,Noverlap,Nfft,r,Wdw)
% [G,Sxdx,Sxx]=STFT_LE(x,delay,Nwind,Noverlap,Nfft,r,Wdw)
% Line Enhancer with short time Fourier transforms (STFT) of length Nfft.
% x is the noisy signal and delay sets its delayed version to be used as 
% a reference for filtering out the predicted output xf.
% G is the estimated frequency filter that is applied on Nwind-long FT of x 
% with Noverlap time samples between them. If specified, the sequences
% are windowed by the 'Wdw' function - default is 'Parzen'.
% r = 1,2,3,... applies a r times shorter window on the short-time sequences x.
% 


if delay < Nwind, disp('Caution : delay should be greater than Nwind for better results'); end

% Windows conditionning
Nwindr = fix(Nwind/r);
if nargin < 7
   Window = parzen(Nwind);
   Windowr = [parzen(Nwindr);zeros(Nwind-Nwindr,1)];
else
   Window = feval(Wdw,Nwind);
   Windowr = feval(Wdw,Nwindr);	Windowr = [Windowr;zeros(Nwind-Nwindr,1)]; 
end
Window = Window/sum(Window);
Windowr = Windowr/sum(Windowr);

% Signals conditioning
xd = x(delay+1:end);	xd = xd(:);
n = length(xd);
x = x(1:n);	x = x(:);

% Ensemble of regressed windows Xd and regressor windows X
Noverlap = ceil(Noverlap);
R = Nwind - Noverlap;		% Shift btw 2 consecutive windows
k = fix((n-Noverlap)/R);	% Number of windows
index = 1:Nwind;
Sxdx = 0;	Sxx = 0;
for i=1:k
   Xdi = fft(Window.*xd(index),Nfft);
   Xi = fft(Windowr.*x(index),Nfft);
   Sxdx = Sxdx + Xdi.*conj(Xi);
   Sxx = Sxx + abs(Xi).^2;
   index = index + R;
end

% Optimal filter 
G = Sxdx./Sxx;










