function Mat = CSM_W(y,x,nfft,Noverlap,Window)     
% Mat = CSM_W(y,x,nfft,Noverlap,Window)
% Welch's estimate of the Cross Spectral Matrix of signals y and x.
%
% Inputs:
% -------
% y, x = matrix of signals, where each column contains a track.
% nfft, Noverlap, and Window are as in function 'PSD' or 'PWELCH' of Matlab.
% Denoting by Nwind the window length, it is recommended to use
% nfft = 2*NWind and Noverlap = 2/3*Nwind with a hanning window or Noverlap = 1/2*Nwind with a half-sine window.
%
% CSM_W calls function 'CPS_W'.
%
% Outputs:
% --------
% S is a structure organized as follows:
%                   Mat.S = Cyclic Cross Spectral Matrix (# columns x)x(# columns y)x(Nwind/2+1)
%                   Mat.f = vector of frequencies
%                   Mat.K = number of blocks
%                   Mat.Var_Reduc = Variance Reduction factor
% 
%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%


Kx = size(x,2);
Ky = size(y,2);

Nw = length(Window);
if Nw == 1
    Nw = Window;
end
Mat.S = zeros(Ky,Kx,Nw/2+1);
count = 1;

same = 0;
if Kx == Ky
    if all(all(y==x))
        same = 1;
    end
end
if same == 1
    for i = 1:Kx
        Sij = CPS_W(y(:,i),x(:,i),nfft,Noverlap,Window);
        Mat.S(i,i,:) = real(Sij.S(1:Nw/2+1));
        for j = 1:i-1
            waitbar(count/Kx/(Kx-1)*2),count = count + 1;
            Sij = CPS_W(y(:,i),x(:,j),nfft,Noverlap,Window);
            Mat.S(i,j,:) = Sij.S(1:Nw/2+1);
            Mat.S(j,i,:) = conj(Sij.S(1:Nw/2+1));
        end
    end
else
    for i = 1:Ky
        for j = 1:Kx
            waitbar(count/Kx/Ky),count = count + 1;
            Sij = CPS_W(y(:,i),x(:,j),nfft,Noverlap,Window);
            Mat.S(i,j,:) = Sij.S(1:Nw/2+1);
        end
    end
end

Mat.Var_Reduc = Sij.Var_Reduc; 
Mat.K = Sij.K;
Mat.f = Sij.f; 


function Spec = CPS_W(y,x,nfft,Noverlap,Window)
% Spec = CPS_W(y,x,Noverlap,Window)
% Welch's estimate of the Cross Power Spectrum of signals y and x.
%
% Inputs:
% -------
% y, x = signals
% nfft, Noverlap, and Window are as in function 'PSD' or 'PWELCH' of Matlab.
% Denoting by Nwind the window length, it is recommended to use
% nfft = 2*NWind and Noverlap = 2/3*Nwind with a hanning window or Noverlap = 1/2*Nwind with a half-sine window.
%
% Outputs
% -------
% Spec is a structure organized as follows:
%                   Spec.S = Cross Power Spectrum vector
%                   Spec.f = vector of normalized frequencies
%                   Spec.K = number of blocks
%                   Spec.Var_Reduc = Variance Reduction factor
%
%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: June 2016
%%%%%%%%%%%%%%%%%%%%%%%%


if length(Window) == 1
    Window = hanning(Window);
end
Window = Window(:);
n = length(x);          % Number of data points
nwind = length(Window); % length of window

% check inputs
if nwind <= Noverlap,error('Window length must be > Noverlap');end
if nfft < nwind,error('Window length must be <= nfft');end

y = y(:);
x = x(:);
K = fix((n-Noverlap)/(nwind-Noverlap));	% Number of windows

% compute CPS
index = 1:nwind;
f = (0:nfft-1)/nfft;
CPS = 0;

for i=1:K
    Yw = fft(Window.*y(index),nfft);		% Yw(f)
    Xw = fft(Window.*x(index),nfft);		% Xw(f-a)
    CPS = Yw.*conj(Xw) + CPS;
    index = index + (nwind - Noverlap);
end

% normalize
% KMU = K*norm(Window)^2;	% Normalizing scale factor ==> asymptotically unbiased
KMU = K;
CPS = CPS/KMU;

% variance reduction factor
Window = Window(:)/norm(Window);
Delta = nwind - Noverlap;
R2w = xcorr(Window);
k = nwind+Delta:Delta:min(2*nwind-1,nwind+Delta*(K-1));
if length(k) >1
    Var_Reduc = R2w(nwind)^2/K + 2/K*(1-(1:length(k))/K)*(R2w(k).^2);
else
    Var_Reduc = R2w(nwind)^2/K;
end

% set up output parameters
Spec.S = CPS;
Spec.f = f;
Spec.K = K;
Spec.Var_Reduc = Var_Reduc;





