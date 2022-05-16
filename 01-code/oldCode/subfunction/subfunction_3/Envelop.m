function a=Envelop(x)
y=abs(hilbert(x));
p=y./sum(y);
p=p+(p==0);
a=-sum(p.*log2(p));