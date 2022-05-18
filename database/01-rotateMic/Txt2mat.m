%TXT2MAT
clc
clear
%close all
subfunction_path1='/Users/wjq/Desktop/ExpDuctMode/Workspace/subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.txt';'*.*'},'mat�����ļ���ȡ','MultiSelect','on');%MultiSelect��ѡ
load([location,'����˵��','/','parameter6000.mat']); %ѡ���ļ���������
if isstr(fname)
   fname=cellstr(fname);
end
for i_file=1:length(fname)

Data = importdata(fullfile(location,char(fname(i_file)))); %ѡ���ļ���������
[Pulse,Rotor_Speed] = keyRotation_RealTime(Data(1:end,end),fs); %ͨ�������źŻ�ȡת����Ϣ.ȥͷȥβ      
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





