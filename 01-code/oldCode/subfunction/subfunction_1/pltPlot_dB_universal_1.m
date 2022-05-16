function [rpm,tsignal,xuhao]=pltPlot_dB_universal(signal,fs,object,objectName,testTime,fname,save_directory,if_RI,type)

[Pulse,Rotor_Speed] = keyRotation_RealTime(signal(1:end,end),fs); %通过键向信号获取转速信息.去头去尾      
H_if_RI=1:(length(Rotor_Speed)-1);

%% 添加失速双锁相功能
%---------------%%-%%%%--------可注释掉-------------%%%%%---------------------
if if_RI ==1
S_data =smooth(signal(:,object(9)),2000,'loess' );%低通滤波%作为失速启动的键向信号
Threshold=mean(S_data);%默认值----->调整参数 
figure
plot(signal(:,object(9)));hold on
plot(S_data);hold on
plot([1 204800],[Threshold Threshold]);
title('确认失速')


disp(['~请再次确认',char(fname),'是否是失速工况？'])
disp(['~需观察是否在频谱图上观察到失速旋转频率和其倍频？如无，非失速工况'])
disp(['~以及确认Threshold值，是否在红线之间？'])
prompt = '观察失速运动情况为1；观察一周变化为0；';
pause

[Pulse_RI,RI_Speed] = RI_RealTime(S_data,fs,Threshold); %通过键向信号获取失速相位信息 

% 每圈都找对应的失速相位
for NumPhase=1:length(Pulse)
    temp_Phase=Pulse_RI./Pulse(NumPhase);
    temp_order=find(temp_Phase<1);
     if isempty(temp_order)
        continue
     end
    if temp_order(end)==length(Pulse_RI)
        break
    end
     
    Phase(NumPhase,:)=(1-(Pulse(NumPhase)-Pulse_RI(temp_order(end)))/(Pulse_RI(temp_order(end)+1)-Pulse_RI(temp_order(end))))*360;   
end
[RI_order(:,1),RI_order(:,2)]=sort(+Phase);%一个符号绝对旋转方向
    H_if_RI=RI_order(:,2)';
end     
%----------------%%%%---------------------%%%%%---------------------
%% 


solutime=1;
for H=H_if_RI  %每转一圈刷新一次网格，一个时间步！%按照RI_marker排序
%% 计算和安排相对位置关系
distance_key_rotor=6;%键向位置和最近叶片之间的距离
R_rotor=184.4002;%mm/动叶半径
x=[2;4.41;7.68;11.14;14.39;17.41;20.73;23.77;26.63;34];
%y_2=[20.29;24.92;29.60;35.37;40.80;44.83;49.71;55.40;59.70;72.36];
y_2=[0;0;0;0;0;0;0;0;0;0];
distance_point_B1_R1=3.8*R_rotor*2*pi/29;
distance_point_B1_C1=5.3/360*R_rotor*2*pi;%17.0575mm;将C1移到70mm，差值为（70-17.0575）mm
distance_pointR1_R1=R_rotor*2*pi/29;
interval=Rotor_Speed(H)/60*R_rotor/fs*2*pi;
round_point=distance_pointR1_R1*29/interval;%一周的点数
Len=floor(round_point)-2;%=130;%调节显示节点个数;取一周围成一个三维的圈
point_Pulse_rotor=round(distance_key_rotor/interval); %%%12000-（19*1.6092）；13000-（12*1.7433）；10000-（5*1.0740）；11000-（-2*1.3598）；125000-（16*1.6763）
point_rotor_rotor=round_point/29; %两个叶片之间的距离来推算采样一个叶道的点数 %144000/BPF
tran_C1=(70-distance_point_B1_C1)/interval;%C1实际点位和理想点位差23.2235个点位，C1前移23个点位,第24
point_B1_R1_interval=round(distance_point_B1_R1./interval);  %B1与R1差52.5761个点位，R1-R8前移53个点位，第54
point_B1_C1_interval=round((y_2(end)-y_2(1)-distance_point_B1_C1)./interval);  %%17.0575mm;将C1移到70mm，差值为（70-17.0575）mm
%% 生成网格矩阵
y_1=(1:Len+1)*interval;
X = repmat(x,size(y_1))';
yy_1 = repmat(y_1,size(x))';
yy_2=repmat(y_2,size(y_1))';
Y=yy_1+yy_2;
% 2维面和3维面转换
initial_y=Y(1,1);
Theta=(Y-initial_y)./R_rotor;
XX=R_rotor*cos(Theta);
YY=R_rotor*sin(Theta);
ZZ=X;

% 导入数据，对网格赋值
    rpm{solutime}(:,1) = signal(Pulse(H):Pulse(H)+Len,object(1));%传感器B1
    rpm{solutime}(:,2) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(2));
    rpm{solutime}(:,3) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(3));
    rpm{solutime}(:,4) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(4));
    rpm{solutime}(:,5) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(5));
    rpm{solutime}(:,6) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(6));
    rpm{solutime}(:,7) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(7));
    rpm{solutime}(:,8) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(8));
    rpm{solutime}(:,9) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(9));
    rpm{solutime}(:,10) = signal((point_B1_C1_interval+Pulse(H)):(point_B1_C1_interval+Pulse(H)+Len),object(10));%传感器C1
   % 一圈有29个叶片，对其进行划分，标记动叶位置（29份）:标记
      for k=1:29  
        Pulse_solid(k)=round(1+point_rotor_rotor*(k-1));%这样就可以在变转速插入叶片了
      end
     %rpm{solutime}(Pulse_solid+point_Pulse_rotor,:)=1.1*max(rpm{solutime}(:,end));  %在大约point_rotor_rotor个点搜寻正确的叶片位置，暂定用方差作为衡量指标

 %%
     tsignal.surfaces(solutime).zonename='mysurface zone';
     tsignal.surfaces(solutime).x=XX;    %size 3x3 
     tsignal.surfaces(solutime).y=ZZ;    %size 3x3
     tsignal.surfaces(solutime).z=YY;    %size 3x3
     tsignal.surfaces(solutime).v(1,:,:)=rpm{solutime};%根据键向信号判读
     tsignal.surfaces(solutime).solutiontime=solutime;
     solutime=solutime+1;
end
 %% 整合数据，生成文件
%title=''; 
rotorSpeed=floor(mean(Rotor_Speed)/100)*100;

if if_RI ==0
NAME = [strrep(char(fname),type,'-'),'-blade1-',date];  %存储文件夹
else
 NAME = ['DoubleLock-',strrep(char(fname),'.mat','-'),'-blade-',date];  %存储文件夹   
end
Varnames=[NAME];
output_file_name=[save_directory,'\',NAME,'-',num2str(rotorSpeed),'.plt']; 
tsignal.Nvar=4;     
tsignal.varnames={'x','y','z',Varnames};
mat2tecplot(tsignal,output_file_name);
%save([save_directory,'\',NAME,'.mat'],'rpm')
%save([save_directory,'\',NAME,'1.mat'],'tsignal');
xuhao=Pulse_solid+point_Pulse_rotor;
%save([save_directory,'\',NAME,'xuhao.mat'],'xuhao');
end
 
 
