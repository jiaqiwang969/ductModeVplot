%Design a minimum-order lowpass filter with a 500 Hz passband cutoff frequency
%and 600 Hz stopband cutoff frequency, with a sampling frequency of 2000 Hz, 
%at least 40 dB attenuation in the stopband, and less than 3 dB of ripple 
%in the passband:

rp = 3;          % Passband ripple
rs = 40;         % Stopband ripple
fs = 2000;       % Sampling frequency
f = [500 600];   % Cutoff frequencies
a = [1 0];       % Desired amplitudes
% Compute deviations
dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)]; 
[n,fo,ao,w] = remezord(f,a,dev,fs);
b = remez(n,fo,ao,w);
freqz(b,1,1024,fs);
title('Lowpass Filter Designed to Specifications');