function frequencyDomainPlot(signal,fs,Scale_data)
    the_freq = [0:length(signal)/Scale_data - 1]*fs/length(signal);  %����Ƶ����ɢ�̶�
    temp_freq = abs(fft(signal))*2/length(signal);
    plot(the_freq, temp_freq(1:length(signal)/Scale_data));%ylabel('A-2(kulite)'); 
   