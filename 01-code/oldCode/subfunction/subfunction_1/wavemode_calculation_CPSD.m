%i_method = 1: 时间谱
%i_method = 2: 阶比谱
function [GAMMA,freq,rotor_speed]=wavemode_calculation_CPSD(data,fs,nk,i_method)
    [Pulse,rotor_speed]=keyRotation(data(:,end),fs);
    
    if i_method == 1
        signal = data(:,1:nk);
    end
    
    if i_method == 2
        Order_Cycle = 2^floor(log2(mean(diff(Pulse))/1.6));   %每转采集点数目
        data_resample = zeros(Order_Cycle*(length(Pulse)-1),nk);
        for j = 1:length(Pulse)-1
            data_resample(Order_Cycle*(j-1)+1:Order_Cycle*j,:) = resample(data(Pulse(j):Pulse(j+1)-1,1:nk), Order_Cycle, Pulse(j+1)-Pulse(j));
        end
        signal = data_resample;
        fs = Order_Cycle;
    end
        
    L_signal = length(signal);
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft = 2^(ceil(log2(L_seg))+1);  
    for k=1:nk
        for l = 1:nk
            [C{k,l},freq] = cpsd(signal(:,k),signal(:,l),Wind,Noverlap,Nfft,fs);          
        end
    end
    GAMMA = zeros(Nfft/2+1,nk+1);
    mode=-nk/2:nk/2;
    for m = 1:nk+1
        temp_f = zeros(Nfft/2+1,1);
        for k = 1:nk
            for l = 1:nk
                temp_f = temp_f + 0.5*C{k,l}*exp(i*mode(m)*2*pi*k/nk)*exp(-i*mode(m)*2*pi*l/nk);
            end
        end
        GAMMA(:,m) = temp_f/(nk*nk);
    end

 end
   
