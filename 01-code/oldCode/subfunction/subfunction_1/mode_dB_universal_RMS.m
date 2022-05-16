function [rotor_speed]=mode_dB_universal_RMS(signal,fs,Scale_data,object,nk,objectName,testTime,fname)
    L_signal = length(signal);
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft = 2^(ceil(log2(L_seg))+1);  
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
    Fk=[1/29,1/2,1,2,3];
    df = fs/Nfft;
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
    %plot(the_freq, imag(temp_freq(:,1:32)));%ylabel('A-2(kulite)'); 

    for k=1:5
        A_mf(k,:)=max(abs(GAMMA(floor(rotor_speed/60*29*Fk(k)/df)+[floor(-2/df):floor(2/df)],:)));
    end
    bar(mode,A_mf');hold on
   legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',16);

   %plot(m,A_mf','.');
    set(gca,'XTick',mode);
    set(gca,'Ygrid','on') 
    title({[testTime,'-模态分析'];[fname,'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
    xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
    %plot(the_freq, imag(temp_freq(:,1:32)));%ylabel('A-2(kulite)'); 
        
 end
   
