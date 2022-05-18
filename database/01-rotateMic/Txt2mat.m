%TXT2MAT
clc
clear
%close all
subfunction_path1='/Users/wjq/Desktop/ExpDuctMode/Workspace/subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.txt';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'参数说明','/','parameter6000.mat']); %选择文件导入数据
if isstr(fname)
   fname=cellstr(fname);
end
for i_file=1:length(fname)

Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
[Pulse,Rotor_Speed] = keyRotation_RealTime(Data(1:end,end),fs); %通过键向信号获取转速信息.去头去尾      
dT=length(Data)/fs;
MeanSpeed=50*round(mean(Rotor_Speed)/50);
% if Rotor_Speed(end)-Rotor_Speed(1)>200
%     label='I';
% else if Rotor_Speed(end)-Rotor_Speed(1)<-200
%          label='D'
% else
%          label='U';
%     end
% end

save([location,'/',strrep(char(fname(i_file)),'.txt',''),'.mat'],'Data');

end





