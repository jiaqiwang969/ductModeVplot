function frequencyDomainPlote(signal,fs)
    the_freq = [0:length(signal)/2.56-1 ]*fs/length(signal);  % ˝æ›∆µ”Ú¿Î…¢øÃ∂»
    temp_freq = abs(fft(signal))*2/length(signal);
    plot(the_freq, temp_freq(1:length(signal)/2.56));%ylabel('A-2(kulite)'); 
    