function a=Envelop2(x,fs)

Nfft = length(x);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_x = hilbert(x);
envelope_spec_e1 = abs(fft(abs(hilbert_x)));
a=max(envelope_spec_e1(find(f>40&f<1000),:));



