%%����������ģ̬��Ƶ�ʽ�ֹ��ϵ��ͼ����������2019-9-26�޸�;
%������FFT�����ļ���Աȣ�����μ�wavemode_calculation_FFT.m����������2019-09-30
%wjq, 2019-10-6�޸ģ��������й���
clc
clear
close all
i_cpsd_method = 2;  %CPSD���������1:ʱ���ף�2���ױ���
i_spectrum_method = 2;    %1:fftδ�ֶ�ƽ����2��fft�ֶ�ƽ����


subfunction_path1='H:\����Ԥ��һ�廯����ͨ�ð�_��ת��ϻ_ver1\subfunction\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat�����ļ���ȡ','MultiSelect','on');%MultiSelect��ѡ
load([location,'\','����˵��','\','parameter.mat']); %ѡ���ļ���������
disp(Note);
% % //======����ͼ����ָ���ļ���===============  
save_directory = [strrep(location,'Database','����������ģ̬��Ƶ�ʽ�ֹ��ϵ��ͼ(CPSD+Spectrum)'),date];  %Ƶ��ͼ�洢�ļ���

%save_directory = ['����������ģ̬��Ƶ�ʽ�ֹ��ϵ��ͼ(CPSD+Spectrum)',date];  %Ƶ��ͼ�洢�ļ���
    
    
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('�ļ��д��ڣ�');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
for i_file=1:length(fname)
    close all
    Data = importdata(fullfile(location,char(fname(i_file)))); %ѡ���ļ���������
    Data=V2Pa_Universal(Data,kulite_transform_ab);
    Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
    
%%  1-Ƶ��ͼ    
    [freq,rotor_speed,freq_dB]=spectrum_calculation_FFT(Data,fs,30,i_spectrum_method,save_directory,fname,i_file);

%% 2-ģ̬ͼ-by CPSD      
    index=[4:15 31];
    [SPECTRA,the_freq]=wavemode_calculation_CPSD(Data(:,index),fs,12,i_cpsd_method,save_directory,fname,i_file);
        
end
   
