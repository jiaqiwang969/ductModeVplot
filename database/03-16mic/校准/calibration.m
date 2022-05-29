clc
clear
[fnameRead,pnameRead] = uigetfile('*.tdms','Select the data file');
fileRead = [pnameRead,fnameRead];
[ConvertedData]=convertTDMS(false,fileRead);
Sound_data = ConvertedData.Data.MeasuredData(3).Data;
Sound_data = Sound_data - mean(Sound_data);
fs = 50000;
Length_data = length(Sound_data);
t = [0:Length_data-1]/fs;
figure
plot(t,Sound_data')
title(fnameRead(1:4))
%%%
data_fft = fs/4;
N_overlap = data_fft/2;
N_seg = round((Length_data-data_fft)/(data_fft-N_overlap)) - 1;
data_freq = zeros(data_fft,1);
for i = 1:N_seg
    data = Sound_data((i-1)*(data_fft-N_overlap)+1:(i-1)*(data_fft-N_overlap)+data_fft,1);
    temp_freq = abs(fft(data))*2/data_fft;
    data_freq = data_freq + temp_freq;
end
data_freq = data_freq/N_seg;
freq_dB = 20*log10(data_freq/2e-5);     % 计算声压级
Scale_data = 5;
data_freq_dB = freq_dB(1:data_fft/Scale_data);
the_freq = [0:data_fft/Scale_data - 1]*fs/data_fft;  %数据频域离散刻度
[maxvalue,indx] = max(freq_dB);
freq_max = the_freq(indx);
figure
plot(the_freq, data_freq_dB); 
xlim([0 5000])
title(fnameRead(1:4))
xlabel(['max\_freq=   ',num2str(freq_max),'  Hz;   ', '  SPL ', num2str(maxvalue),' dB'])
SPL = 94;  %标准声源 dB
Sensitivity_ini = 1;  %每通道的初始灵敏度，mv/pa;
Standard_pa = 10^(SPL/20)*(2e-5);
Actual_pa = 10^(maxvalue/20)*(2e-5);
Sensitivity_actual = Actual_pa*Sensitivity_ini/Standard_pa


