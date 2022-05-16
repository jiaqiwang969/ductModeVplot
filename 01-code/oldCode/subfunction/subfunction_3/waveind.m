function s=waveind(x)
%to caculate the wave index of the signal
s=rootmeansquare(x)./meanamp(x);