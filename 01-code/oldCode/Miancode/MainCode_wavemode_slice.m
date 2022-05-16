%选择特定频率下的环形声阵列模态输出通用程序，董广明，2019-9-27
%增加了FFT方法的计算对比，具体参见wavemode_calculation_FFT.m，董广明，2019-09-30
clc
clear
close all
% frequency selection
i_sync = 1;  %选择同步频率，否则置零
i_cpsd = 1;   %利用互功率谱计算，否则置零
i_fft_method = 3;    %FFT计算参数：1:fft未分段平均；2：fft分段平均；3：fft阶比谱
i_cpsd_method = 1;    %cpsd计算参数：1:fft未分段平均；2：fft阶比谱
i_order_spec = 1;  %进行阶比分析，否则置零
if i_order_spec ==1
    i_fft_method = 3;
end

if i_sync == 1
    filename_label = 'Sync';
    Freq_slice = [1/29,1/2,1,2,3];    %对应BPF
else
    Freq_slice = [2320, 2450];       %对应Hz
    filename_label = [num2str(Freq_slice), ' Hz'];
end

if i_cpsd == 0
    filename_label =[filename_label, ' --FFT method ', num2str(i_fft_method)];
end

subfunction_path1='F:\动画预测一体化程序通用版_ver3\subfunction\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'\','参数说明','\','parameter.mat']); %选择文件导入数据
disp(Note);
% % //======保存图像至指定文件夹===============  

if i_cpsd ==1
    save_directory = [location,'\','特定频率下的环形阵列声模态(CPSD)',date];  %频谱图存储文件夹
else
    save_directory = [location,'\','特定频率下的环形阵列声模态(FFT)',date];  %频谱图存储文件夹
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
mode = [-32/2:32/2];
for i_file=1:length(fname)
    Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
    Data=V2Pa_Universal(Data,kulite_transform_ab);
    Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
    if i_cpsd ==1 
        [SPECTRA,freq,rotor_speed]=wavemode_calculation_CPSD(Data,fs,32,i_cpsd_method);
    else
        [SPECTRA,freq,rotor_speed]=wavemode_calculation_FFT(Data,fs,32,i_fft_method);
    end
    df =freq(2) - freq(1);
    h=figure('Visible', 'on');
    set(gcf,'outerposition',get(0,'screensize'));%最大化
    
    if i_sync == 1                  %转频同步频率
        if i_order_spec == 1     %阶次
            for k=1:length(Freq_slice)
                Wavemode(k,:)=max(abs(SPECTRA(floor(29*Freq_slice(k)/df)+[floor(-0.1/df):floor(0.1/df)],:)));
            end
        else                              %频率
            for k=1:length(Freq_slice)
                Wavemode(k,:)=max(abs(SPECTRA(floor(rotor_speed/60*29*Freq_slice(k)/df)+[floor(-5/df):floor(5/df)],:)));
            end
        end
        bar(mode,Wavemode');hold on
        legend({'1*SSF';'1/2*BPF';'1*BPF';'2*BPF';'3*BPF';},'Location','NorthEast','FontSize',12);
    else
        if i_order_spec == 1
            for k=1:length(Freq_slice)
                Wavemode(k,:)=max(abs(SPECTRA(floor(Freq_slice(k)*60/(df*rotor_speed))+[floor(-0.1/df):floor(0.1/df)],:)));
            end
        else
            for k=1:length(Freq_slice)
                Wavemode(k,:)=max(abs(SPECTRA(floor(Freq_slice(k)/df)+[floor(-5/df):floor(5/df)],:)));
            end
        end
        bar(mode,Wavemode');hold on
        legend(cellstr(num2str(Freq_slice')),'Location','NorthEast','FontSize',12);       
    end
    set(gca,'XTick',mode);
    set(gca,'Ygrid','on') 
    if i_cpsd == 1
        title({[testTime,'-模态分析'];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
    else
        title({[testTime,'-模态分析', ' -FFT method ', num2str(i_fft_method)];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
    end
    xlabel('Mode Number：m','FontSize',16);ylabel('Amplitude','FontSize',16);
    
    saveas(h,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.fig'])
    saveas(h,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.png'])
    close all
end

