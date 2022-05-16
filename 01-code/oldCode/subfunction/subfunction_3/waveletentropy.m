function a=waveletentropy(X)
X=abs(X);
p=X./sum(X);
p=p+(p==0);
a=-sum(p.*log2(p))./log2(length(X));