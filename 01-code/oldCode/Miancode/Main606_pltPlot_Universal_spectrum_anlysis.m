
clc
clear
close all
subfunction_path1='D:\研究\research\Compressor Sound\2019-压气机动态信号测试\试验17-2019-11-12-ROTATE\subfunction\subfunction_1';
addpath(subfunction_path1);
fntype='.mat';
[fname,location]=uigetfile({['*',fntype];'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'\','参数说明','\','parameter.mat']); %选择文件导入数据
disp(Note);
% % //======保存图像至指定文件夹===============  
save_directory = [strrep(location,'Database','Frequency_DMD')]; 
%save_directory = [strrep(location,'Database','Frequency_DMD_DATA12')];  %频谱图存储文件夹
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
Data(:,1:end-1)=Data(:,1:end-1);%-mean(Data(:,1:end-1));
[rpm,tsignal,xuhao,Rotor_Speed]=pltPlot_dB_universal(Data,fs,object,objectName,testTime,char(fname(i_file)),save_directory,0,fntype);
[the_freq,freq_dB]=frequencyDomainPlot_dB(Data,fs,2.56);
close all
%[the_freq,freq_dB]=frequencyDomainPlot_dB_no_deal(Data,fs,2.56);   
h=figure;set(gcf,'outerposition',get(0,'screensize'));%最大化
[tsignal2,Info]=computeDMD(rpm,tsignal,save_directory,char(fname(i_file)),mean(Rotor_Speed),fs,the_freq,freq_dB,fntype);

h1 = figure('Color',[1 1 1]);set(gcf,'outerposition',get(0,'screensize'));axes6 = axes('Parent',h1)
object_1=[object];
    for k=1:length(object_1)
    hs1 = plot3(the_freq, (k*0.7).*ones(size(the_freq)), freq_dB(:,object_1(k)),'Parent',axes6,'Color',[0 0 0]);
    hold on
    end

xlabel(axes6,{'Norm. Frequency (f/f_r_o_t)'});
ylabel(axes6,{'Sensor Array'});
zlabel(axes6,'Power Spectrum (dB)');
xlim(axes6,[15 mean(Rotor_Speed)/60*29*1.2]);
ylim(axes6,[0.69,7]);
zlim(axes6,[105,160]);
view(axes6,[0.100000000000273 81.9999999999999]);
box(axes6,'on');
set(axes6,'XGrid','on'); 
set(axes6,'FontSize',10,'XGrid','on','YTick',...
    [0.7 1.4 2.1 2.8 3.5 4.2 4.9 5.6 6.3 7],'YTickLabel',...
    {'B1','R1','R2','R3','R4','R5','R6','R7','R8','C1'});
bname=[];
label=unique(Info.DMD_mode(:,2));
for k=1:length(unique(Info.DMD_mode(:,2)))
aname{1,k}=num2str(label(k));
end
set(axes6,'FontSize',10,'XGrid','on','XTick',[10*round(mean(Rotor_Speed)/60*label/10)]',...
    'XTickLabel',aname)
for k=1:length(Info.DMD_mode(:,7))
hold on
plot3([Info.DMD_mode(k,7),Info.DMD_mode(k,7),Info.DMD_mode(k,7),Info.DMD_mode(k,7),Info.DMD_mode(k,7),Info.DMD_mode(k,7),Info.DMD_mode(k,7),Info.DMD_mode(k,7)],...
    [0.69,0.69,0.69,7,7,7,0.69,7],[105,160,105,105,105,160,160,160],'Parent',axes6);
end
%arrow([2450 0.1 126.4],[2450 0.4 126.4],25,'BaseAngle',60)
title({[strrep(char(fname(i_file)),fntype,'-'),num2str(round(mean(Rotor_Speed))),'rpm-',num2str(fs),'-Av10Hz'];['First-six-mode:',num2str(Info.DMD_mode(1:6,2)')]}  )
saveas(h1,[save_directory,'\','h1',strrep(char(fname(i_file)),fntype,'-'),...
    num2str(round(mean(Rotor_Speed))),'rpm-',num2str(fs),'.fig'])
saveas(h1,[save_directory,'\','h1',strrep(char(fname(i_file)),fntype,'-'),...
     num2str(round(mean(Rotor_Speed))),'rpm-',num2str(fs),'.png'])
% pause
% close all
 
end


