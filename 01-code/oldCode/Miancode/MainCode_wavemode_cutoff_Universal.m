%%环形声阵列模态与频率截止关系绘图，董广明，2019-9-26修改
%增加了FFT方法的计算对比，具体参见wavemode_calculation_FFT.m，董广明，2019-09-30

clc
clear
close all

i_order_spec = 1;  %进行阶比分析，否则置零
i_cpsd = 1;   %利用互功率谱计算，否则置零
i_fft_method = 3;    %1:fft未分段平均；2：fft分段平均；3：fft阶比谱
i_cpsd_method = 2;  %CPSD计算参数：1:时间谱；2：阶比谱
if i_order_spec ==1
    i_fft_method = 3;
    i_cpsd_method = 2;
end

% subfunction_path1='D:\动画预测一体化程序通用版_ver3\subfunction\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'\','参数说明','\','parameter.mat']); %选择文件导入数据
disp(Note);
% % //======保存图像至指定文件夹===============  
if i_cpsd ==1
    save_directory = [location,'\','环形声阵列模态与频率截止关系绘图(CPSD)',date];  %频谱图存储文件夹
    filename_label = ['CPSD method ', num2str(i_cpsd_method)];
else
    save_directory = [location,'\','环形声阵列模态与频率截止关系绘图(FFT)',date];  %频谱图存储文件夹
    filename_label =['FFT method ', num2str(i_fft_method)];
end
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
for i_file=1:length(fname)
    Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
    Data=V2Pa_Universal(Data,kulite_transform_ab);
    Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
    if i_cpsd ==1 
        [SPECTRA,freq,rotor_speed]=wavemode_calculation_CPSD(Data,fs,32, i_cpsd_method);
        SPECTRA = 10*log10(abs(SPECTRA)/4e-10);
%         SPECTRA = abs(SPECTRA);
    else
        [SPECTRA,freq,rotor_speed,data_freq]=wavemode_calculation_FFT(Data,fs,32,i_fft_method);
        SPECTRA=20*log10(abs(SPECTRA)/2e-5); 
    end
    h=figure('Visible', 'on');
    set(gcf,'outerposition',get(0,'screensize'));%最大化
    imagesc([-32/2:32/2],freq,SPECTRA); 
    axis xy

    if i_order_spec == 1
        ylim([1,30]);
        if i_cpsd == 1
            title({[testTime,'-截止模态分析',' -CPSD method ', num2str(i_cpsd_method)];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
        else
            title({[testTime,'-截止模态分析', ' -FFT method ', num2str(i_fft_method)];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
        end
        xlabel('Mode Number：m','FontSize',16);ylabel('Order','FontSize',16);        
    else
         ylim([1,rotor_speed*2]);
        if i_cpsd == 1
            title({[testTime,'-截止模态分析',' -CPSD method ', num2str(i_cpsd_method)];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
        else
            title({[testTime,'-截止模态分析', ' -FFT method ', num2str(i_fft_method)];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
        end
        xlabel('Mode Number：m','FontSize',16);ylabel('Frequency (Hz)','FontSize',16);
    end
    df=freq(2)-freq(1);
    [x1,y1]=find(SPECTRA==max(max(SPECTRA(round(10/df):round(20/df),:))),1);Zdata=SPECTRA(x1,:);Zdata(y1)=0;[y1_1]=find(Zdata==max(Zdata));
    list=[round(0.5/df) round(1.5/df);round(14/df) round(14.5/df);round(28.5/df) round(29.5/df);round(58.5/df) round(59.5/df);round(10/df) round(20/df);round(10/df) x1-20;x1+20 round(20/df)];
    listName={'SSF';'12BPF';'1BPF';'2BPF';'RI1';'RI2';'RI3';};
    for k=1:7 
    [x(k),y(k)]=find(SPECTRA==max(max(SPECTRA(list(k,1):list(k,2),:))),1);Zdata=SPECTRA(x(k),:);Zdata(y(k))=0;[y_1(k)]=find(Zdata==max(Zdata),1);
    o_RI(k)=x(k)*df;m_RI(k)=y(k)-17;m_RI_1(k)=y_1(k)-17;
    text(m_RI(k),o_RI(k),{[num2str(round(m_RI(k)*10)/10),',',num2str(round(SPECTRA(x(k),y(k))))];[num2str(o_RI(k)),',',num2str(round(o_RI(k)*rotor_speed/60))]});text(m_RI_1(k),o_RI(k),[num2str(m_RI_1(k)),',',num2str(round(SPECTRA(x(k),y_1(k))))]);   
    end
    saveas(h,[save_directory,'\','Image','-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.fig'])
    saveas(h,[save_directory,'\','Image','-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.png'])

    for k=1:7
    h0{k}=figure;bar([-32/2:32/2],SPECTRA(x(k),:));hold on
    ylim([70 125]);
    title({[listName{k},'-',char(fname(i_file))];['f=',num2str(rotor_speed/60*o_RI(k)),',',num2str(round(o_RI(k)*rotor_speed/60)),';Order=',num2str(o_RI(k)),',m=',num2str(m_RI(k))]},'FontSize',14)
    saveas(h0{k},[save_directory,'\',num2str(k),'-',listName{k},'-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.fig'])
    saveas(h0{k},[save_directory,'\',num2str(k),'-',listName{k},'-',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.png'])
    end
    close all
end
   
