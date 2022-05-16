%evaluation of theoretical frequency gain of ALE filter
%evaluation of 3dB bandwidth and side-lobe relative amplitude
clc
clear
for iwind = 1:1:15
    Nwind = 2^(iwind+3);
    Nfft = 16*Nwind;
    windowname = 'rectwin';
    w = feval(windowname,Nwind);
    w = w/sum(w);
    W = fft(w,Nfft);
    SNR = 1;
    H_abs = (0.5*SNR*Nwind*W.*conj(W))./(1+0.5*SNR*Nwind*W.*conj(W));
    H_abs_db = convert2db(H_abs(1:Nfft/2));
    f = [0:Nfft/2-1]/Nfft;
    MLWidth_1(iwind) = 2*f(max(find((H_abs_db-H_abs_db(1))>=-3)));
    FSLHeight_1(iwind) = max(H_abs(round(Nfft/Nwind):Nfft/2))/H_abs(1);
end
figure
plot(2.^[4:1:18],MLWidth_1)
hold on

for iwind = 1:1:15
    Nwind = 2^(iwind+3);
    Nfft = 16*Nwind;
    windowname = 'hanning';
    w = feval(windowname,Nwind);
    w = w/sum(w);
    W = fft(w,Nfft);
    SNR = 1;
    H_abs = (0.5*SNR*Nwind*W.*conj(W))./(1+0.5*SNR*Nwind*W.*conj(W));
    H_abs_db = convert2db(H_abs(1:Nfft/2));
    f = [0:Nfft/2-1]/Nfft;
    MLWidth_2(iwind) = 2*f(max(find((H_abs_db-H_abs_db(1))>=-3)));
    FSLHeight_2(iwind) = max(H_abs(round(2*Nfft/Nwind):Nfft/2))/H_abs(1);
end
plot(2.^[4:1:18],MLWidth_2,'r')
hold on

for iwind = 1:1:15
    Nwind = 2^(iwind+3);
    Nfft = 16*Nwind;
    windowname = 'parzenwin';
    w = feval(windowname,Nwind);
    w = w/sum(w);
    W = fft(w,Nfft);
    SNR = 1;
    H_abs = (0.5*SNR*Nwind*W.*conj(W))./(1+0.5*SNR*Nwind*W.*conj(W));
    H_abs_db = convert2db(H_abs(1:Nfft/2));
    f = [0:Nfft/2-1]/Nfft;
    MLWidth_3(iwind) = 2*f(max(find((H_abs_db-H_abs_db(1))>=-3)));
    FSLHeight_3(iwind) = max(H_abs(round(4*Nfft/Nwind):Nfft/2))/H_abs(1);
end
plot(2.^[4:1:18],MLWidth_3,'g')
hold on

for iwind = 1:1:15
    Nwind = 2^(iwind+3);
    SNR = 1;
    Nfft = 16*Nwind;
    f = [0:Nfft/2-1]/Nfft;
    Diric_abs = abs(diric(2*pi*f,Nwind));
    H_abs = (0.5*SNR*Nwind*Diric_abs)/(1+0.5*SNR*Nwind);   
    H_abs_db = convert2db(H_abs(1:Nfft/2));
    MLWidth(iwind) = 2*f(max(find((H_abs_db-H_abs_db(1))>=-3)));
    FSLHeight_4(iwind) = max(H_abs(round(Nfft/Nwind):Nfft/2))/H_abs(1);
end
plot(2.^[4:1:18],MLWidth,'k')

figure
plot(2.^[4:1:18],MLWidth_1./MLWidth)
hold on
plot(2.^[4:1:18],MLWidth_2./MLWidth,'r')
hold on
plot(2.^[4:1:18],MLWidth_3./MLWidth,'g')

figure
plot(2.^[4:1:18],FSLHeight_1)
hold on
plot(2.^[4:1:18],FSLHeight_2,'r')
hold on
plot(2.^[4:1:18],FSLHeight_3,'g')
hold on
plot(2.^[4:1:18],FSLHeight_4,'k')
