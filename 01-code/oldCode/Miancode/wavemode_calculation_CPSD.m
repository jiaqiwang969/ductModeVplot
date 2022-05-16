%i_method = 1: 时间谱
%i_method = 2: 阶比谱
function [GAMMA,freq]=wavemode_calculation_CPSD(data,fs,nk,i_method,save_directory,fname,i_file)
    [Pulse,rotor_speed]=keyRotation(data(:,end),fs);
    filename_label = ['CPSD method ', num2str(i_method)];
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
    
    GAMMA = 10*log10(abs(GAMMA)/4e-10);
%       GAMMA = abs(GAMMA);
    h=figure('Visible', 'on');
    set(gcf,'outerposition',get(0,'screensize'));%最大化
    imagesc([-nk/2:nk/2],freq,GAMMA); 
    axis xy; ylim([1,30]);
    testTime='试验17-2019-11-11';
    title({[testTime,'-截止模态分析',' -CPSD method ', num2str(i_method)];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
    xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);
    
    df=freq(2)-freq(1);
    [x1,y1]=find(GAMMA==max(max(GAMMA(round(10/df):round(20/df),:))),1);Zdata=GAMMA(x1,:);Zdata(y1)=0;[y1_1]=find(Zdata==max(Zdata));
    list=[round(0.5/df) round(1.5/df);round(14/df) round(14.5/df);round(28.5/df) round(29.5/df);round(58.5/df) round(59.5/df);round(10/df) round(20/df);round(10/df) x1-20;x1 round(20/df)];
    listName={'SSF';'12BPF';'1BPF';'2BPF';'RI1';'RI2';'RI3';};
    for k=1:7 
    [x(k),y(k)]=find(GAMMA==max(max(GAMMA(list(k,1):list(k,2),:))),1);Zdata=GAMMA(x(k),:);Zdata(y(k))=0;[y_1(k)]=find(Zdata==max(Zdata),1);
    o_RI(k)=x(k)*df;m_RI(k)=y(k)-(nk/2+1);m_RI_1(k)=y_1(k)-(nk/2+1);
    text(m_RI(k),o_RI(k),{[num2str(round(m_RI(k)*10)/10),',',num2str(round(GAMMA(x(k),y(k))))];[num2str(o_RI(k)),',',num2str(round(o_RI(k)*rotor_speed/60))]});text(m_RI_1(k),o_RI(k),[num2str(m_RI_1(k)),',',num2str(round(GAMMA(x(k),y_1(k))))]);   
    end
    saveas(h,[save_directory,'\','Image','-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.fig'])
    saveas(h,[save_directory,'\','Image','-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.png'])

    for k=1:7
    h1{k}=figure;bar([-nk/2:nk/2],GAMMA(x(k),:));hold on
    ylim([70 150]);
    title({[listName{k},'-',char(fname(i_file))];['f=',num2str(rotor_speed/60*o_RI(k)),',',num2str(round(o_RI(k)*rotor_speed/60)),';Order=',num2str(o_RI(k)),',m=',num2str(m_RI(k))]},'FontSize',14)
    saveas(h1{k},[save_directory,'\',num2str(k),'-',listName{k},'-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.fig'])
    saveas(h1{k},[save_directory,'\',num2str(k),'-',listName{k},'-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.png'])
    end

 end
   
