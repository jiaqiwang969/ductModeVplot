%ѡ���ض�Ƶ���µĻ���������ģ̬���ͨ�ó��򣬶�������2019-9-27
%������FFT�����ļ���Աȣ�����μ�wavemode_calculation_FFT.m����������2019-09-30
clc
clear
close all
% frequency selection
i_sync = 1;  %ѡ��ͬ��Ƶ�ʣ���������
i_cpsd = 1;   %���û������׼��㣬��������
i_fft_method = 3;    %FFT���������1:fftδ�ֶ�ƽ����2��fft�ֶ�ƽ����3��fft�ױ���
i_cpsd_method = 1;    %cpsd���������1:fftδ�ֶ�ƽ����2��fft�ױ���
i_order_spec = 1;  %���нױȷ�������������
if i_order_spec ==1
    i_fft_method = 3;
end

if i_sync == 1
    filename_label = 'Sync';
    Freq_slice = [1/29,1/2,1,2,3];    %��ӦBPF
else
    Freq_slice = [2320, 2450];       %��ӦHz
    filename_label = [num2str(Freq_slice), ' Hz'];
end

if i_cpsd == 0
    filename_label =[filename_label, ' --FFT method ', num2str(i_fft_method)];
end

subfunction_path1='F:\����Ԥ��һ�廯����ͨ�ð�_ver3\subfunction\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat�����ļ���ȡ','MultiSelect','on');%MultiSelect��ѡ
load([location,'\','����˵��','\','parameter.mat']); %ѡ���ļ���������
disp(Note);
% % //======����ͼ����ָ���ļ���===============  

if i_cpsd ==1
    save_directory = [location,'\','�ض�Ƶ���µĻ���������ģ̬(CPSD)',date];  %Ƶ��ͼ�洢�ļ���
else
    save_directory = [location,'\','�ض�Ƶ���µĻ���������ģ̬(FFT)',date];  %Ƶ��ͼ�洢�ļ���
end
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('�ļ��д��ڣ�');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
mode = [-32/2:32/2];
for i_file=1:length(fname)
    Data = importdata(fullfile(location,char(fname(i_file)))); %ѡ���ļ���������
    Data=V2Pa_Universal(Data,kulite_transform_ab);
    Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
    if i_cpsd ==1 
        [SPECTRA,freq,rotor_speed]=wavemode_calculation_CPSD(Data,fs,32,i_cpsd_method);
    else
        [SPECTRA,freq,rotor_speed]=wavemode_calculation_FFT(Data,fs,32,i_fft_method);
    end
    df =freq(2) - freq(1);
    h=figure('Visible', 'on');
    set(gcf,'outerposition',get(0,'screensize'));%���
    
    if i_sync == 1                  %תƵͬ��Ƶ��
        if i_order_spec == 1     %�״�
            for k=1:length(Freq_slice)
                Wavemode(k,:)=max(abs(SPECTRA(floor(29*Freq_slice(k)/df)+[floor(-0.1/df):floor(0.1/df)],:)));
            end
        else                              %Ƶ��
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
        title({[testTime,'-ģ̬����'];[char(fname(i_file)),'-ת��: ',num2str(rotor_speed),'-�����ʣ�',num2str(fs)]},'FontSize',14)
    else
        title({[testTime,'-ģ̬����', ' -FFT method ', num2str(i_fft_method)];[char(fname(i_file)),'-ת��: ',num2str(rotor_speed),'-�����ʣ�',num2str(fs)]},'FontSize',14)
    end
    xlabel('Mode Number��m','FontSize',16);ylabel('Amplitude','FontSize',16);
    
    saveas(h,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.fig'])
    saveas(h,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-',filename_label,'.png'])
    close all
end

