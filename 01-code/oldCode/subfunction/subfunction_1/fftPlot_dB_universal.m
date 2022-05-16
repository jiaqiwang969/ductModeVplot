function [rotor_speed,the_freq,freq_dB_scale]=fftPlot_dB_universal(signal,fs,Scale_data,object,objectName,testTime,fname)
    figure
    data_fft = fs/10;
    N_overlap = data_fft/2;
    N_seg = round((length(signal)-data_fft)/(data_fft-N_overlap)) - 1;
    data_freq = zeros(data_fft,size(signal,2));
    for i = 1:N_seg
        data = signal((i-1)*(data_fft-N_overlap)+1:(i-1)*(data_fft-N_overlap)+data_fft,:);
        temp_freq = abs(fft(data))*2/data_fft;
        data_freq = data_freq + temp_freq;
    end
    data_freq = data_freq/N_seg;
    freq_dB =20*log10(data_freq./2e-5);     % 计算声压级
    the_freq = [0:data_fft/Scale_data - 1]*fs/data_fft;  %数据频域离散刻度
    freq_dB=freq_dB(1:data_fft/Scale_data,:);
    freq_dB_scale=freq_dB(1:data_fft/Scale_data,:);
    plot(the_freq, freq_dB_scale(:,:)); hold on
    legend(objectName,'Location','SouthEast');
   
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
    ik = floor(fs/(Scale_data*rotor_speed/60*29));
    for nk=1:ik
      plot(nk*rotor_speed/60*29,min(min(freq_dB))-1,'^');
      hold on
    end
     set(gca,'XLim',[0,3.5*rotor_speed/60*29+1])
     title({[testTime,'-FFT频谱分析'];[fname,'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
     xlabel('频率/Hz','FontSize',16);ylabel('幅值/dB','FontSize',16);
 end
   
