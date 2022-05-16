
clc
clear
%close all
subfunction_path1='D:\研究\research\Compressor Sound\CODE_wjq\subfunction\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.txt';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'\','参数说明','\','parameter.mat']); %选择文件导入数据
disp(Note);
% % //======保存图像至指定文件夹===============  
save_directory = [location,'\','46个传感器的频谱图',date];  %频谱图存储文件夹
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
Data(:,1:end-1)=Data(:,1:end-1);%-mean(Data(:,1:end-1));

h1=figure
[the_freq,freq_dB]=frequencyDomainPlot_dB(Data,fs,2.56);
[Pulse,rotor_speed]=keyRotation(Data(:,end),fs);
set(gcf,'color',[1 1 1],'position',[10 10 770 520]);
    
h2 =figure('Visible', 'on');
set(gcf,'outerposition',get(0,'screensize'));%最大化
axes1 = axes('Parent',h2);
hold(axes1,'on');
object_1=[1:32];
    for k=object_1
    hs1 = plot3(the_freq, (k*0.4).*ones(size(the_freq)), freq_dB(:,k)+40,'-k');
    hold on
    end
object_2=[1:14];
    for k=object_2
    hs1 = plot3(the_freq, 32*0.4+(k*0.7).*ones(size(the_freq)), freq_dB(:,32+k),'-k');
    hold on
    end

xlabel({'Norm. Frequency (f/f_r_o_t)'});
ylabel({'Sensor Array'});
zlabel('Power Spectrum (dB)');
xlim(axes1,[15 6000]);
ylim(axes1,[0.69,32*0.4+14*0.7]);
zlim(axes1,[105,155]);
view(axes1,[9.70000000000021 82.7999999999999]);
box(axes1,'on');
set(axes1,'XGrid','on');   
set(axes1,'FontSize',14,'XGrid','on','XTick',[mean(rotor_speed)/60*[1:29]],...
    'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29'},'YTick',...
    [0.4*[1:32] 32*0.4+([1:14]*0.7)],'YTickLabel',...
    {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','A1','A2','A3','A4','B1','B2','B3','B4','C1','C2','C3','C4','D1','D2'});
title({[testTime,'-频谱分析', ' -平均-32(+40dB) '];[char(fname(i_file)),'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
saveas(h1,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-','1D','.fig'])
saveas(h1,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-','1D','.png'])
saveas(h2,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-','2D','.fig'])
saveas(h2,[save_directory,'\',strrep(char(fname(i_file)),'.mat','-'),num2str(rotor_speed),'rpm-','2D','.png'])
end

