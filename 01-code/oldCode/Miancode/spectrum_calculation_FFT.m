%wjq,2019-10-06
%自功率谱@abs then average
% i_method = 1, direct FFT without average; 
% i_method = 2, FFT with data segmentation and average; 

function [the_freq,rotor_speed,freq_dB]=spectrum_calculation_FFT(signal,fs,nk,i_method,save_directory,fname,i_file)
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
     filename_label = ['FFT ', num2str(i_method)];

    if i_method == 1
        data_fft = 2^floor(log2(length(signal)));
        data = signal(1:data_fft,1:nk);
        the_freq = [0:data_fft/2.56 - 1]*fs/data_fft;  %数据频域离散刻度
        data_freq = fft(data)*2/data_fft;
        data_freq = data_freq(1:data_fft/2.56,:); 
    end
    if i_method == 2
        data_fft = fs/10;
        N_overlap = data_fft/2;
        N_seg = round((length(signal)-data_fft)/(data_fft-N_overlap)) - 1;
        data_freq = zeros(data_fft,1);
        for i_seg = 1:N_seg
            data = signal((i_seg-1)*(data_fft-N_overlap)+1:(i_seg-1)*(data_fft-N_overlap)+data_fft,1:nk);
            temp_freq = abs(fft(data))*2/data_fft;
            data_freq = data_freq + temp_freq;
        end
        data_freq=data_freq(1:data_fft/2.56,:)/N_seg;
        the_freq = [0:data_fft/2.56 - 1]*fs/data_fft;  %数据频域离散刻度
    end
    freq_dB=20*log10(abs(data_freq)/2e-5);
    
    h0=figure('Visible', 'on');
    set(gcf,'outerposition',get(0,'screensize'));%最大化
    for k=1:3
        plot3(the_freq, ((k)*0.7).*ones(size(the_freq)), freq_dB(:,(k))+60,'-g');hold on       
    end  
    for k=4:16
        plot3(the_freq, ((k)*0.7).*ones(size(the_freq)), freq_dB(:,(k)),'-r');hold on       
    end
    for k=17:19
        plot3(the_freq, ((k)*0.7).*ones(size(the_freq)), freq_dB(:,(k)),'-k');hold on       
    end
    for k=21:30
        plot3(the_freq, ((k)*0.7).*ones(size(the_freq)), freq_dB(:,(k)),'-b');hold on       
    end 
    for k=20
        plot3(the_freq, ((k)*0.7).*ones(size(the_freq)), freq_dB(:,(k)),'-b');hold on       
    end 
    
     xlabel({'Norm. Frequency (f/f_r_o_t)'});ylabel({'Sensor Array'});zlabel('Power Spectrum (dB)');
     xlim([15 1.1*rotor_speed/60*29]);ylim([0.69,nk*0.7]);zlim([70,160]);
     view([9.70000000000021 82.7999999999999]);
     testTime='试验17-2019-11-11';
     title({[testTime,'-自功率谱分析',' -FFT method ', num2str(i_method)];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
    
    saveas(h0,[save_directory,'\','Spectrum-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.fig'])
    saveas(h0,[save_directory,'\','Spectrum-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.png'])
  
 end
   
