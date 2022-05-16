%%环形声阵列模态与频率截止关系绘图，董广明，2019-9-26修改;
%增加了FFT方法的计算对比，具体参见wavemode_calculation_FFT.m，董广明，2019-09-30
%wjq, 2019-10-6修改，集成所有功能
clc
clear
close all
i_cpsd_method = 2;  %CPSD计算参数：1:时间谱；2：阶比谱
i_spectrum_method = 2;    %1:fft未分段平均；2：fft分段平均；


subfunction_path1='H:\动画预测一体化程序通用版_旋转机匣_ver1\subfunction\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'\','参数说明','\','parameter.mat']); %选择文件导入数据
disp(Note);
% % //======保存图像至指定文件夹===============  
save_directory = [strrep(location,'Database','环形声阵列模态与频率截止关系绘图(CPSD+Spectrum)'),date];  %频谱图存储文件夹

%save_directory = ['环形声阵列模态与频率截止关系绘图(CPSD+Spectrum)',date];  %频谱图存储文件夹
    
    
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
    close all
    Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
    Data=V2Pa_Universal(Data,kulite_transform_ab);
    Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
    
%%  1-频谱图    
    [freq,rotor_speed,freq_dB]=spectrum_calculation_FFT(Data,fs,30,i_spectrum_method,save_directory,fname,i_file);

%% 2-模态图-by CPSD      
    index=[4:15 31];
    [SPECTRA,the_freq]=wavemode_calculation_CPSD(Data(:,index),fs,12,i_cpsd_method,save_directory,fname,i_file);
        
end
   
