%add comparison 2019-09-30 dgm
% i_method = 1, direct FFT without average; 
% i_method = 2, FFT with data segmentation and average; 
% i_method = 3, FFT after constant angle resampling; 

function [a_mf,the_freq,rotor_speed]=wavemode_calculation_FFT(signal,fs,nk,i_method)
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
    
    if i_method == 1
        data_fft = 2^floor(log2(length(signal)));
        data = signal(1:data_fft,1:32);
        the_freq = [0:data_fft/2.56 - 1]*fs/data_fft;  %����Ƶ����ɢ�̶�
        data_freq = fft(data)*2/data_fft;
        data_freq = data_freq(1:data_fft/2.56,:); 
    end
    
    if i_method == 2
        data_fft = 2^floor(log2(fs/5));
        the_freq = [0:data_fft/2.56 - 1]*fs/data_fft;  %����Ƶ����ɢ�̶�
        N_overlap = data_fft/2;
        N_seg = round((length(signal)-data_fft)/(data_fft-N_overlap)) - 1;
        data_freq = zeros(data_fft,32);
        for i_seg = 1:N_seg
            data = signal((i_seg-1)*(data_fft-N_overlap)+1:(i_seg-1)*(data_fft-N_overlap)+data_fft,1:32);
            temp_freq = fft(data)*2/data_fft;
            data_freq = data_freq + temp_freq;
        end
        data_freq=data_freq(1:data_fft/2.56,:)/N_seg;   
    end
    
    if i_method == 3
        data = signal(:,1:32);
        Order_Cycle = 2^floor(log2(mean(diff(Pulse))/1.6));   %ÿת�ɼ�����Ŀ
        data_resample = zeros(Order_Cycle*(length(Pulse)-1),32);
        for j = 1:length(Pulse)-1
            data_resample(Order_Cycle*(j-1)+1:Order_Cycle*j,:) = resample(data(Pulse(j):Pulse(j+1)-1,:), Order_Cycle, Pulse(j+1)-Pulse(j));
        end
        data_fft = length(data_resample); %�ȽǶ��ز���������ݳ���
        the_freq = [0:data_fft/2.56-1]*Order_Cycle/data_fft;
        data_freq = fft(data_resample)*2/data_fft;
        data_freq = data_freq(1:data_fft/2.56,:); 
    end
    
    m=-nk/2:nk/2;
    for k=1:length(m)
        a_mf(:,k)=1/nk*data_freq(:,1:nk)*exp(m(k)*i*2*pi*(1:nk)/nk).'; 
    end
 end
   
