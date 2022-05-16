function deconstruction_signal_1(signal,fs,object)
fk=12000;
x=resample(signal(:,1),fk,fs);%qian

fn=104;
%band-pass filtration

%f=[2*0.85*fn/fk 2*0.9*fn/fk 2*1.1*fn/fk 2*1.15*fn/fk];% Cutoff frequencies
f=[2*10/fk 2*20/fk 2*300/fk 2*320/fk];
A=[0 1 0];% Desired amplitudes
rp=0.153;% Passband ripple
rs=16.92;% Stopband ripple
devp=1-10^(-rp/20);
devs=10^(-rs/20);
dev=[devp devs devp];
[n,f0,A0,w]=remezord(f,A,dev);

if rem(n,2)
   n=n+1;
end

b=remez(n,f0,A0,w);

[h,w]=freqz(b,1,length(b),1);
hr=abs(h);
h=abs(h);
h=20*log10(h);
figure
subplot(211)
plot(b);grid;
subplot(212)
plot(w,h);grid;

x_filter = filter(b,1,x);  

%STFT_LE

delay = 10;
Nwind = 1024;
NFFT = 2*Nwind;
N = length(x_filter);

% version with absolute magnitude TF
G = STFT_LE(x_filter,Nwind+delay,Nwind,0,NFFT,1,'hanning');
G1 = abs(G);
[y1,g] = Filt_STFT(x_filter,G1.');    %复数转置G1.'; 复数共轭转置G1';
N1 = length(y1);
e1 = x_filter(1:N1,:) - y1;
% version with complex TF
G2 = G;
y2 = Filt_STFT(x_filter,G2.');
N2 = Nwind + delay;
e2 = x_filter(N2+1:N) - y2(1:N-N2);
    
t = (1:N);
figure,plot(x_filter),hold on,plot(y1,'r'),plot(t(N2+1:N),y2(1:N-N2),':g')
xlabel('y')
figure,plot(e1,'r'),hold on,plot(t(N2+1:N),e2(1:N-N2),':g')
xlabel('e')
[psd_x,f_x] = pwelch(x_filter,[],[],[],fs);
[psd_e1,f_e1] = pwelch(e1,[],[],[],fs);
[psd_e2,f_e2] = pwelch(e2,[],[],[],fs);
figure
plot(f_x,psd_x);
hold on
plot(f_e1,psd_e1,'r');
plot(f_e2,psd_e2,':g');
xlabel('x and e')
[psd_y1,f_y1] = pwelch(y1,[],[],[],fs);
[psd_y2,f_y2] = pwelch(y2,[],[],[],fs);
figure
plot(f_x,psd_x);
hold on
plot(f_y1,psd_y1,'r');
plot(f_y2,psd_y2,':g');
xlabel('x and y')

% Demodulation

Nfft = length(x);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_x = hilbert(x);
envelope_spec_x = abs(fft(abs(hilbert_x)));
figure
plot(f,envelope_spec_x(1:Nfft/2)./max(envelope_spec_x(3:end)))
hold on
Nfft = length(x_filter);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_x_filter = hilbert(x_filter);
envelope_spec_x_filter = abs(fft(abs(hilbert_x_filter)));
plot(f,envelope_spec_x_filter(1:Nfft/2)./max(envelope_spec_x_filter(3:end)),'k')
hold on
Nfft = length(e1);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_e1 = hilbert(e1);
envelope_spec_e1 = abs(fft(abs(hilbert_e1)));
plot(f,envelope_spec_e1(1:Nfft/2)./max(envelope_spec_e1(3:end)),'r')
hold on
Nfft = length(e2);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_e2 = hilbert(e2);
envelope_spec_e2 = abs(fft(abs(hilbert_e2)));
plot(f,envelope_spec_e2(1:Nfft/2)./max(envelope_spec_e2(3:end)),'g')


end
 
      


   
