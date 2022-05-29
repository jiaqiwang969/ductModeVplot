function [rotor_speed]=fftPlot_dB_universal(signal,fs,Scale_data,object,objectName,testTime,fname)

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
    plot(the_freq, freq_dB(1:data_fft/Scale_data,object)); hold on
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
     
     


Nstd = 0.2;%噪声标准差
NR = 500;%采集分次
MaxIter = 5000;%最大迭代次数
lag=1;

[modes its]=ceemdan(signal(:,object(1)),Nstd,NR,MaxIter,lag);
t=1:length(ecg);

[a b]=size(modes);

figure;
subplot(a+1,1,1);
plot(t,ecg);% the ECG signal is in the first row of the subplot
ylabel('ECG')
set(gca,'xtick',[])
axis tight;

for i=2:a
    subplot(a+1,1,i);
    plot(t,modes(i-1,:));
    ylabel (['IMF ' num2str(i-1)]);
    set(gca,'xtick',[])
    xlim([1 length(ecg)])
end;

subplot(a+1,1,a+1)
plot(t,modes(a,:))
ylabel(['IMF ' num2str(a)])
xlim([1 length(ecg)])

figure;
boxplot(its);
     
 end
   
