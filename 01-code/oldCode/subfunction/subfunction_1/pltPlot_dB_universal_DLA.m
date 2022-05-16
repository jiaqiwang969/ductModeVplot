function [RPM1,tsignal,Pulse_blade,T,signal_X]=pltPlot_dB_universal_ver1(signal,fs,object,objectName,testTime,fname,save_directory,if_RI,type,average_degree)

[Pulse,Rotor_Speed] = keyRotation_RealTime(signal(1:end,end),fs); %通过键向信号获取转速信息.去头去尾      
N_blade=4;%添加一个功能-3或4个叶道为一组-增大数据组数
for kkk=1:floor(29/N_blade)-1
Pulse(1:end-1,kkk+1)=Pulse(1:end-1,1)+round(kkk*(diff(Pulse(:,1))/29*N_blade));
Rotor_Speed(1:end-1,kkk+1)=Rotor_Speed(1:end-1,1);
end
Pulse_temp=Pulse(1:end-1,:)';
Rotor_Speed_temp=Rotor_Speed(1:end-1,:)';
Pulse_divide=Pulse_temp(:);
Rotor_Speed_divide=Rotor_Speed_temp(:);
Len=round(min(diff(Pulse(:,1))));%/29*N_blade


H_if_RI=1:(length(Rotor_Speed_divide)-1);

%% 添加失速双锁相功能
%---------------%%-%%%%--------可注释掉-------------%%%%%---------------------
if if_RI ==1
S_data =smooth(signal(:,object(4)),1500,'loess' );%低通滤波%作为失速启动的键向信号
Threshold=mean(S_data);%默认值----->调整参数 
figure('InvertHardcopy','off','Color',[1 1 1]);

axes1 = axes;
hold(axes1,'on');
plot1=plot(signal(:,object(4)));
set(plot1(1),'Color',[0.862745106220245 0.862745106220245 0.862745106220245]);
hold on
plot(S_data,'LineWidth',3);hold on
plot([1 204800],[Threshold Threshold],'LineStyle','--');
title('确认失速')


disp(['~请再次确认',char(fname),'是否是失速工况？'])
disp(['~需观察是否在频谱图上观察到失速旋转频率和其倍频？如无，非失速工况'])
disp(['~以及确认Threshold值，是否在红线之间？'])
prompt = '观察失速运动情况为1；观察一周变化为0；';
[Pulse_RI,Pulse_RI_low1,Pulse_RI_high1,RI_Speed] = RI_PhaseLockPoint(S_data,fs,Threshold); %锁相点
hold on
plot(Pulse_RI,S_data(Pulse_RI),'MarkerSize',25,'Marker','.','LineWidth',15,'LineStyle','none')
hold on
plot(Pulse_RI_low1,S_data(Pulse_RI_low1),'MarkerSize',25,'Marker','.','LineStyle','none','Color',[1 0 0])
xlim(axes1,[Pulse_RI(2) Pulse_RI(6)]);
plot(signal(:,object(4)))
%pause

T=1/(mean(Rotor_Speed_divide)-mean(RI_Speed))*60;

% 每圈都找对应的失速相位
for NumPhase=1:length(Pulse_divide)
    temp_Phase=Pulse_RI_low1./Pulse_divide(NumPhase);
    temp_order=find(temp_Phase<1);
     if isempty(temp_order)
        continue
     end
    if temp_order(end)==length(Pulse_RI_low1)
        break
    end
     
    Phase_divide(NumPhase,:)=(1-(Pulse_divide(NumPhase)-Pulse_RI_low1(temp_order(end)))/(Pulse_RI_low1(temp_order(end)+1)-Pulse_RI_low1(temp_order(end))))*360;   
end
    [RI_order(:,1),RI_order(:,2)]=sort(+Phase_divide);%一个符号绝对旋转方向
    H_if_RI=RI_order(:,2)';

end 


%----------------%%%%---------------------%%%%%---------------------
%% 
%添加双锁相“平均”模块

solutime=1;
for H=H_if_RI  %每转一圈刷新一次网格，一个时间步！%按照RI_marker排序
%% 计算和安排相对位置关系

distance_key_rotor=6;%键向位置和最近叶片之间的距离
R_rotor=184.4002;%mm/动叶半径
x=[2;4.41;7.68;11.14;14.39;17.41;20.73;23.77;26.63;34];
y_2=[20.29;24.92;29.60;35.37;40.80;44.83;49.71;55.40;59.70;72.36];
distance_point_B1_R1=3.8*R_rotor*2*pi/29;
distance_point_B1_C1=5.3/360*R_rotor*2*pi;%17.0575mm;将C1移到70mm，差值为（70-17.0575）mm
distance_pointR1_R1=R_rotor*2*pi/29;
interval=Rotor_Speed_divide(H)/60*R_rotor/fs*2*pi;
round_point=distance_pointR1_R1*29/interval;%一周的点数
%Len=floor(round_point)-2;%=130;%调节显示节点个数;取一周围成一个三维的圈
point_Pulse_rotor=round(distance_key_rotor/interval); %%%12000-（19*1.6092）；13000-（12*1.7433）；10000-（5*1.0740）；11000-（-2*1.3598）；125000-（16*1.6763）
point_rotor_rotor=round_point/29; %两个叶片之间的距离来推算采样一个叶道的点数 %144000/BPF
tran_C1=(70-distance_point_B1_C1)/interval;%C1实际点位和理想点位差23.2235个点位，C1前移23个点位,第24
point_B1_R1_interval=round(distance_point_B1_R1./interval);  %B1与R1差52.5761个点位，R1-R8前移53个点位，第54
point_B1_C1_interval=round((y_2(end)-y_2(1)-distance_point_B1_C1)./interval);  %%17.0575mm;将C1移到70mm，差值为（70-17.0575）mm


% 导入数据，对网格赋值
    rpm{solutime}(:,1) = signal(Pulse_divide(H):Pulse_divide(H)+Len,object(1));%传感器B1
    rpm{solutime}(:,2) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(2));
    rpm{solutime}(:,3) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(3));
    rpm{solutime}(:,4) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(4));
    rpm{solutime}(:,5) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(5));
    rpm{solutime}(:,6) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(6));
    rpm{solutime}(:,7) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(7));
    rpm{solutime}(:,8) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(8));
    rpm{solutime}(:,9) = signal((point_B1_R1_interval+Pulse_divide(H)):(point_B1_R1_interval+Pulse_divide(H)+Len),object(9));
    rpm{solutime}(:,10)= signal((point_B1_C1_interval+Pulse_divide(H)):(point_B1_C1_interval+Pulse_divide(H)+Len),object(10));%传感器C1
  
    solutime=solutime+1;

end
 %% 整合数据，生成文件
%title=''; 
rotorSpeed=floor(mean(Rotor_Speed_divide)/100)*100;
signal_X(:,1)=signal(1:end-point_B1_R1_interval+1,object(1));
signal_X(:,2)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,3)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,4)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,5)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,6)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,7)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,8)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,9)=signal(point_B1_R1_interval:end,object(2));
signal_X(:,10)=signal(point_B1_C1_interval:end-point_B1_R1_interval+point_B1_C1_interval,object(2));

for k=1:360/average_degree
    %% 生成网格矩阵
y_1=(1:Len+1)*interval;
X = repmat(x,size(y_1))';
yy_1 = repmat(y_1,size(x))';
yy_2=repmat(y_2,size(y_1))';
Y=yy_1+yy_2;
Z=zeros(size(X));
% % 2维面和3维面转换
initial_y=Y(1,1);
Theta=(Y-initial_y)./R_rotor;
XX=R_rotor*cos(Theta);
YY=R_rotor*sin(Theta);
ZZ=X;
    
    Xuhao{k}=find(RI_order(:,1)>average_degree*(k-1)&RI_order(:,1)<average_degree*(k));
    rpm_temp=zeros(Len+1,10);
    for kk=1:length(Xuhao{k})
        rpm_temp=rpm_temp+rpm{Xuhao{k}(kk)};
    end
    RPM{k}=rpm_temp/kk;  
    RPM1{k}=RPM{k};
    
     % 一圈有29个叶片，对其进行划分，标记动叶位置（29份）:标记
      for kkk=1:29  
        Pulse_solid(kkk)=round(1+point_rotor_rotor*(kkk-1));%这样就可以在变转速插入叶片了
      end
     Pulse_blade= Pulse_solid+point_Pulse_rotor;
     %RPM{k}(Pulse_blade,:)=1.1*max(RPM{k}(:,end));  %在大约point_rotor_rotor个点搜寻正确的叶片位置，暂定用方差作为衡量指标

 %%
     tsignal.surfaces(k).zonename='mysurface zone';
     tsignal.surfaces(k).x=XX;    %size 3x3 
     tsignal.surfaces(k).y=ZZ;    %size 3x3
     tsignal.surfaces(k).z=YY;    %size 3x3
%     tsignal.surfaces(k).x=Y;    %size 3x3 
%     tsignal.surfaces(k).y=X;    %size 3x3
%     tsignal.surfaces(k).z=Z;    %size 3x3
     tsignal.surfaces(k).v(1,:,:)=RPM{k};%根据键向信号判读
     tsignal.surfaces(k).solutiontime=k;
     
end

if if_RI ==0
NAME = [strrep(char(fname),type,'-'),'-blade-d',num2str(average_degree)];  %存储文件夹
else
NAME = ['DoubleLock-',strrep(char(fname),'.mat','-'),'-blade-d',num2str(average_degree)];  %存储文件夹   
end
Varnames=[NAME];
output_file_name=[save_directory,'\',NAME,'-',num2str(rotorSpeed),'.plt']; 
tsignal.Nvar=4;     
% tsignal.varnames={'x','y','z',Varnames};
tsignal.varnames={'x','y','z',Varnames};
mat2tecplot(tsignal,output_file_name);
%save([save_directory,'\',NAME,'.mat'],'rpm')
%save([save_directory,'\',NAME,'1.mat'],'tsignal');
%save([save_directory,'\',NAME,'xuhao.mat'],'xuhao');
end
 
      


   
