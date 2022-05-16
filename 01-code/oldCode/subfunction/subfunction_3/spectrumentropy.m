function a=spectrumentropy(x)
X=abs(fft(x));
p=X./sum(X);
p=p+(p==0);
a=-sum(p.*log2(p))./log2(length(x));