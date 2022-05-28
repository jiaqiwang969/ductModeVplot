function [S,alpha,f,STFT,t,Nv] = Fast_SC(x,Nw,alpha_max,Fs,opt)
% [S,alpha,f] = Fast_SC(x,Nw,alpha_max,Fs,opt)
% Fast estimation of the Spectral Correlation based on the Short-Time-Fourier-Transform (STFT).
% (stand-alone Matlab code)
%
% Inputs:
% -------
% x: signal whose spectral correlation is to be computed.
% Nw: window size used in STFT.
% alpha_max: upper bound of the cyclic frequency range to scan (in Hz).
% Fs: sampling frequency in Hz (default = 1Hz).
% opt.coh = 1: the spectral COHERENCE is computed instead of the (default) spectral CORRELATION.
%
% Outputs:
% --------
% S: spectral correlation matrix ; dimension = (positive frequencies)x(cyclic frequencies).
% alpha: vector of cyclic frequencies in Hz.
% f: vector of spectral frequencies in Hz.
% STFT: Short-Time-Fourier-Transform matrix ; dimension = (positive frequencies)x(time instants).
% t: vector of time instants of the STFT (in seconds).
% Nv: number of overlapped samples.
%
% Recommendations for setting the paramaters:
% -------------------------------------------
% 1) Nw sets the frequency resolution df ~ 1.5*Fs/Nw (in Hz). The shorter Nw,
% the faster the algorithm. A power of 2 is recommended so as to make
% possible the use of the FFT.
% 2) alpha_max should be set to the largest cyclic frequency of interest,
% but not larger. The shorter it is, the faster the algorithm.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: J. Antoni, September 2016
% Reference: 
% J. Antoni, G. Xin, N. Hamzaoui, "Fast Computation of the Spectral Correlation", 
% Mechanical Systems and Signal Processing, Elsevier, 2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check inputs
if alpha_max > Fs/2
    error('''alpha_max'' must be smaller than Fs/2!');
end
if alpha_max < 0
    error('''alpha_max'' must be non-negative!');
end

% Set value of overlap
% ====================
[Nv,dt] = param_Fast_SC(length(x),Nw,alpha_max,Fs);

% Computation of short-time Fourier transform 
% ===========================================
[STFT,f,t] = LiteSpectrogram(x,Nw,Nv,Nw,Fs);

% Fast spectral correlation/coherence
% ===================================
[S,alpha] = Fast_SC_STFT(STFT,dt,Nw,Fs,[],opt);

I = find(alpha <= alpha_max);
alpha = alpha(I);
S = S(:,I);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nv,dt,da,df] = param_Fast_SC(L,Nw,alpha_max,Fs)
% [Nv,dt,da,df] = param_Fast_SC(L,Nw,alpha_max,Fs)
% Computes the number of overlapped samples Nv to be used in the Fast STFT-based
% estimator of the Spectral Correlation given a signal length L, a window length Nw, 
% the maximum cyclic frequency alpha_max, and the sampling frequency Fs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% block shift
R = fix(Fs/2/alpha_max);	
R = max(1,min(R,fix(.25*Nw)));

% block overlap
Nv = Nw - R;  

% time resolution of STFT (in s)
dt = R/Fs;   

if nargout > 2
    % cyclic frequency resolution (in Hz)
    da = Fs/L;
    
    if nargout > 3
        % carrier frequency resolution (in Hz)
        df = Fs/Nw*sum(Nw.^2)/mean(Nw)^2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,alpha,W,Winv] = Fast_SC_STFT(STFT,Dt,Wind,Fs,Nfft,opt)
% [S,alpha,W] = Fast_SC(STFT,Dt,Window,Fs,Nfft)
% Fast estimation of the Spectral Correlation by scanning several zooms obtained from the short-time-Fourier-transform.
%
% Inputs:
% -------
% STFT: Short-Time-Fourier-Transform matrix ; dimension = (positive frequencies)x(time instants)
% as obtained from a standard algorithm, such as 'spectrogram.m'.
% (if STFT has 3 dimensions, the CROSS-spectral correlation is computed between STFT(:,:,1) and STFT(:,:,2))
% Dt: sampling period (in seconds) of the STFT.
% Window: window used in STFT.
% Fs: sampling frequency in Hz (default = 1Hz).
% Nfft: number of FFT bins used for cyclic frequencies (should be greater than size(STFT,2)).
% opt.abs = 0 or 1: summation of the complex values (default) or absolute magnitudes of the scans.
% opt.calib = 1 or 0: S is equalized by the spectral "scanning" window (default) or not.
% opt.coh = 1: the spectral COHERENCE is computed instead of the (default) spectral CORRELATION.
%
% Outputs:
% --------
% S: "scanned" spectral correlation matrix ; dimension = (positive frequencies)x(cyclic frequencies).
% alpha: vector of cyclic frequencies in Hz.
% W: spectral "scanning window".
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: August 2015
% Revised: July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%

NF = size(STFT,1);
Nw = 2*(NF-1);      % window length

flag = 0;
if nargin < 4
    Fs = 1;
end
if nargin < 5
    Nfft = size(STFT,2);
elseif isempty(Nfft)
    Nfft = size(STFT,2);
end
if nargin < 6
    opt.abs = 0;
    opt.calib = 1;
    opt.coh = 0;
else
    if ~isfield(opt,'abs')
        opt.abs = 0;
    end
    if ~isfield(opt,'calib')
        opt.calib = 1;
    end
    if ~isfield(opt,'coh')
        opt.coh = 0;
    end
end

if length(Wind) == 1
    Wind = hanning(Wind);
end

% -------------------------------------------------------------------------
% Whitening the STFT for computing the spectral coherence
if opt.coh == 1
    Sx = mean(abs(STFT).^2,2);  % Mean power spectral density
    STFT = STFT.*repmat(1./sqrt(Sx),1,size(STFT,2));
end

% -------------------------------------------------------------------------
% Computation of the cyclic modulation spectrum
[S,alpha] = CPS_STFT_zoom(0,STFT,Dt,Wind,Fs,Nfft,flag);
W0 = Window_STFT_zoom(alpha,0,Dt,Wind,Nfft,Fs,'full');
% % S = abs(S);
% % W = abs(W0);
if opt.abs == 1
    S = abs(S);
    W = abs(W0);
else
    W = W0;
end
W(ceil(Nfft/2)+1:Nfft) = 0; % truncate negative frequencies

% -------------------------------------------------------------------------
% Number of scans
Fa = 1/Dt;               % cyclic sampling frequency in Hz
K = fix(Nw/2*Fa/Fs);
% K = Nw/2-1;   % maximum number of scans!!

% -------------------------------------------------------------------------
% Computation and summation of "cross-frequency" cyclic modulation spectra
for k = 1:K
    % positive cyclic frequencies
    [Stemp,alpha,alpha0] = CPS_STFT_zoom(k/Nw*Fs,STFT,Dt,Wind,Fs,Nfft,flag);
    Wtemp = Shift_Window_STFT_zoom(W0,alpha0/Fa*Nfft,'trunc');
    % Wtemp = Window_STFT_zoom(alpha,alpha0,Dt,Wind,Nfft,Fs,'trunc');
    if opt.abs == 1
        S(:,2:Nfft) = S(:,2:Nfft) + abs(Stemp(:,2:Nfft));
        W(2:Nfft) = W(2:Nfft) + abs(Wtemp(2:Nfft));
    else
        S(:,2:Nfft) = S(:,2:Nfft) + Stemp(:,2:Nfft);
        W(2:Nfft) = W(2:Nfft) + Wtemp(2:Nfft);
    end
    %plot(abs(Wtemp),':'),plot(W)
    
    % negative cyclic frequencies
    [Stemp,alpha,alpha0] = CPS_STFT_zoom(-k/Nw*Fs,STFT,Dt,Wind,Fs,Nfft,flag);
    Wtemp = Shift_Window_STFT_zoom(W0,alpha0/Fa*Nfft,'trunc');
    % Wtemp = Window_STFT_zoom(alpha,alpha0,Dt,Wind,Nfft,Fs,'trunc');
    if opt.abs == 1
        S(:,2:Nfft) = S(:,2:Nfft) + abs(Stemp(:,2:Nfft));
        W(2:Nfft) = W(2:Nfft) + abs(Wtemp(2:Nfft));
    else
        S(:,2:Nfft) = S(:,2:Nfft) + Stemp(:,2:Nfft);
        W(2:Nfft) = W(2:Nfft) + Wtemp(2:Nfft);
    end
    %plot(abs(Wtemp),':'),plot(W)
end

% -------------------------------------------------------------------------
% Calibration
if opt.calib == 1
    % Winv = W./(W.^2 + size(STFT,2)/size(STFT,1));   % regularized inverse to avoid division by zero
    % (size(STFT,2)/size(STFT,1) is a rough evaluation of the estimation noise in S)
    % NSR = Fs/(Nw*size(STFT,2)*Dt);
    % Winv = W./(W.^2 + 1e-2*max(abs(W))*sqrt(NSR));
    % S = S.*repmat(Winv.',NF,1)*sum(Wind(:).^2);
    Winv = ones(Nfft,1);
    I = find(W < .5*W(1));
    Winv(1:I(1)) = 1./W(1:I(1));
    Winv(I(1)+1:Nfft) = 1/W(I(1)+1);
    S = S.*repmat(Winv.',NF,1)*sum(Wind(:).^2);
else
    Winv = 1/W(1);
    S = S*Winv*sum(Wind(:).^2);
end
% Impose real values at zero cyclic frequency
% S(:,1) = real(S(:,1));

if opt.coh == 1
    S = S/mean(S(:,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines of Fast_SC_STFT.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,alpha,alpha0,fk,Fa] = CPS_STFT_zoom(alpha0,STFT,Dt,Window,Fs,Nfft,flag)
% [S,alpha,W,alpha0,fk,Fa] = CPS_STFT_zoom(alpha0,STFT,Dt,Window,Fs,Nfft,flag)
% Fast estimation of the zoomed Spectral Correlation from the Short-Time-Fourier-Transform.
%
% Inputs:
% -------
% alpha0: zoomed cyclic frequency (in Hz).
% STFT: Short-Time-Fourier-Transform matrix (positive frequencies)x(times) as obtained from a standard algorithm
% (if STFT has 3 dimensions, the cross-spectral correlation is computed between STFT(:,:,1) and STFT(:,:,2)).
% Dt: sampling period (in seconds) of the STFT.
% Window: window used in STFT.
% Fs: sampling frequency in Hz (default = 1Hz).
% Nfft: number of FFT bins used for cyclic frequencies (should be greater than size(STFT,2)).
%
% Outputs:
% --------
% S: zoomed Spectral Correlation matrix ; dimension = (positive frequencies)x(cyclic frequencies).
% alpha: vector of cyclic frequencies in Hz.
% W: spectral "zooming" window.
% alpha0: actual zoomed cyclic frequency in Hz (may be different from input value due to round-off to an integer number of frequency bins).
% fk: frequency bin number corresponding to alpha0.
% Fa: cyclic sampling frequency (in Hz) (Fa/2 = maximum cyclic frequency).
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: August 2015
%%%%%%%%%%%%%%%%%%%%%%%%%

[NF,NT,N3] = size(STFT);
Nw = 2*(NF-1);          % window length

Fa = 1/Dt;    % cyclic sampling frequency in Hz

if nargin < 5
    Fs = 1;
end
if nargin < 6
    Nfft = NT;
elseif Nfft < NT
    error('Nfft must be greater than or equal to the number of time samples in STFT!')
end

% -------------------------------------------------------------------------
% Check for aliasing
if nargin < 7
    if abs(alpha0) > Fa/2
        disp(['|alpha0| must be selected smaller than ',num2str(Fa/2),'!!']);
    end
end

% -------------------------------------------------------------------------
% Vector of cyclic frequencies
alpha = (0:Nfft-1)/Nfft*Fa;

% -------------------------------------------------------------------------
% Computation "cross-frequency" cyclic modulation spectrum
fk = round(alpha0/Fs*Nw);
alpha0 = fk/Nw*Fs;
if N3 == 1
    if fk >= 0
        S = [STFT(1+fk:NF,:);zeros(fk,NT)].*conj(STFT);
    else
        S = [conj(STFT(1-fk:NF,:));zeros(-fk,NT)].*STFT;
    end
else
    if fk >= 0
        S = [squeeze(STFT(1+fk:NF,:,1));zeros(fk,NT)].*conj(squeeze(STFT(:,:,2)));
    else
        S = [conj(squeeze(STFT(1-fk:NF,:,1)));zeros(-fk,NT)].*squeeze(STFT(:,:,2));
    end
end
S = fft(S,Nfft,2)/NT;

% -------------------------------------------------------------------------
% Calibration
S = S/sum(Window(:).^2)/Fs;

% -------------------------------------------------------------------------
% Removal of aliased cyclic frequencies
ak = round(alpha0/Fa*Nfft);
S(:,ceil(Nfft/2)+1+ak:Nfft) = 0;

% -------------------------------------------------------------------------
% Phase correction
[~,Iw] = max(Window);
% S = S.*repmat(exp(-2i*pi*Iw*(alpha-alpha0-alpha(ceil(NT/2)))/Fs),NF,1);
S = S.*repmat(exp(-2i*pi*Iw(1)*(alpha-alpha0)/Fs),NF,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,W1,W2] = Window_STFT_zoom(alpha,alpha0,Dt,Window,Nfft,Fs,opt)
% W = Window_STFT_zoom(alpha,alpha0,Dt,Window,Nfft,Fs,opt)
% Returns the spectral "zooming window" of the STFT-based zoom algorithm.
%
% Inputs:
% -------
% alpha: vector of cyclic frequencies in Hz.
% alpha0: zoomed cyclic frequency in Hz.
% Dt: sampling period (in seconds) of the STFT.
% Window: window used in STFT.
% Fs: sampling frequency in Hz (default = 1Hz).
% Nfft: number of points in FFT.
% opt: the full window is returned ('full') or its truncated version ('trunc') without aliasing.
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: August 2015
%%%%%%%%%%%%%%%%%%%%%%%%%

Fa = 1/Dt;    % cyclic sampling frequency in Hz

if nargin < 6
    Fs = 1;
end

% -------------------------------------------------------------------------
% Computation the "zooming" window
WSquared = Window(:).^2;
[~,Iw] = max(Window); % set origin of time to the centre of symmetry (maximum value) of the window
W1 = zeros(Nfft,1);
W2 = zeros(Nfft,1);
n = (1:Iw-1)'/Fs;
for k = 1:Nfft
    % "positive" frequencies
    W1(k) = WSquared(Iw) + 2*sum(WSquared(Iw-1:-1:1).*cos(2*pi*n*(alpha(k)-alpha0)));
    % "negative" frequencies (aliased)
    W2(k) = WSquared(Iw) + 2*sum(WSquared(Iw-1:-1:1).*cos(2*pi*n*(alpha(k)-alpha0-Fa)));
end
W = W1 + W2;
% Note: sum(W2) = max(W)

% -------------------------------------------------------------------------
% Removal aliased cyclic frequencies
if strcmp(opt,'trunc')
    ak = round(alpha0/Fa*Nfft);
    W(ceil(Nfft/2)+1+ak:Nfft) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = Shift_Window_STFT_zoom(W0,a0,opt)
% W = Shift_Window_STFT_zoom(W0,a0,opt)
% Circular shift of window W0 by a0 (possibly non-integer) samples.
% This function uses fast linear interpolation.
% opt: the full window is returned ('full') or its truncated version ('trunc') without aliasing.
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: August 2015
%%%%%%%%%%%%%%%%%%%%%%%%%

Nfft = length(W0);

% -------------------------------------------------------------------------
% Circular shift with linear interpolation for non-integer shifts
a1 = floor(a0);
a2 = ceil(a0);
if a1 == a2
    W = circshift(W0,a0);
else
    W = circshift(W0,a1)*(1-(a0-a1)) + circshift(W0,a2)*(a0-a1);
end

% -------------------------------------------------------------------------
% Removal of aliased cyclic frequencies
if strcmp(opt,'trunc')
    W(ceil(Nfft/2)+1+round(a0):Nfft) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [X,f,t] = LiteSpectrogram(x,Window,Noverlap,Nfft,Fs)     
% [X,f,t] = LiteSpectrogram(x,Window,Noverlap,Nfft)
% Computation of the STFT of signal x.
% This is a substantially alleviated and faster version of the 'SPECTROGRAM' code of Matlab.
% Input arguments {Window, Noverlap,Nfft,Fs} and output arguments {X,f,}
% are as in function 'SPECTROGRAM' of Matlab2015.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni: November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Window) == 1
    Window = hanning(Window);
end
Window = Window(:);
n = length(x);          % Number of data points
nwind = length(Window); % length of window
R = nwind - Noverlap;   % block shift

x = x(:);		
K = fix((n-Noverlap)/(nwind-Noverlap));	% Number of windows

% compute STFT
index(1) = 1;
index(2) = nwind;

X = zeros(Nfft/2+1,K);
for k = 1:K
    Xw = fft(Window.*x(index(1):index(2)),Nfft);		% Xw(f-a/2) or Xw(f-a)
    X(:,k) = Xw(1:Nfft/2+1);
    index = index + R;
end

if nargout > 1
    if nargin < 5
        Fs = 1;
    end
    f = (0:Nfft/2)/Nfft*Fs;
    t = (nwind/2:R:nwind/2+(K-1)*R)/Fs;
end








