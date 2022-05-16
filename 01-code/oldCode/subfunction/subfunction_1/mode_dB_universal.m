function [rotor_speed]=mode_dB_universal(signal,fs,Scale_data,object,nk,objectName,testTime,fname)
    the_freq = [0:length(signal)/2.56 - 1]*fs/length(signal);  %数据频域离散刻度
    temp_freq =fft(signal(:,:))*2./length(signal);
    temp_freq=temp_freq(1:length(signal)/2.56,:);
    m=-30:30;
    for k=1:length(m)
    a_mf(:,k)=1/nk*temp_freq(:,1:nk)*exp(m(k)*i*2*pi*(1:nk)/nk).'; 
    %aaa=exp(2*i*2*pi(1:32)/32);
    %bbb=exp(2*i*2*pi(1:32)/32)';
    %ccc=exp(2*i*2*pi(1:32)/32).';
    %a_mf_dB=20*log10(abs(a_mf)./2e-5);
    end
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
    Fk=[1/29,1/2,1,2,3];
    for k=1:5
        A_mf(k,:)=max(abs(a_mf(floor(rotor_speed/60*29*Fk(k)/(fs/length(signal)))+[floor(-2/(fs/length(signal))):floor(2/(fs/length(signal)))],:)));
    end
    bar(m,A_mf');hold on
   legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',16);

   %plot(m,A_mf','.');
    set(gca,'XTick',m);
    set(gca,'Ygrid','on') 
    title({[testTime,'-模态分析'];[fname,'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
    xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
    %plot(the_freq, imag(temp_freq(:,1:32)));%ylabel('A-2(kulite)'); 
    
        
 end
   
