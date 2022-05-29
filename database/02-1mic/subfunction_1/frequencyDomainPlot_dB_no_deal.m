function [the_freq,freq_dB]=frequencyDomainPlot_dB_no_deal(signal,fs,Scale_data)

    the_freq = [0:length(signal)/Scale_data - 1]*fs/length(signal);  %数据频域离散刻度
    data_freq = abs(fft(signal))*2/length(signal);
    freq_dB =20*log10(data_freq/2e-5); 
    hold on% 计算声压级
    freq_dB=freq_dB(1:length(signal)/Scale_data,:);
    h1=plot(repmat(the_freq',1,size(signal,2)), freq_dB);%ylabel('A-2(kulite)'); 

end