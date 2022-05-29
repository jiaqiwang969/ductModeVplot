function [rpm,tsignal,xuhao]=pltPlot_dB_universal(signal,fs,object,objectName,testTime,fname,save_directory,if_RI,type)

[Pulse,Rotor_Speed] = keyRotation_RealTime(signal(1:end,end),fs); %ͨ�������źŻ�ȡת����Ϣ.ȥͷȥβ      
H_if_RI=1:(length(Rotor_Speed)-1);

%% ���ʧ��˫���๦��
%---------------%%-%%%%--------��ע�͵�-------------%%%%%---------------------
if if_RI ==1
S_data =smooth(signal(:,object(9)),2000,'loess' );%��ͨ�˲�%��Ϊʧ�������ļ����ź�
Threshold=mean(S_data);%Ĭ��ֵ----->�������� 
figure
plot(signal(:,object(9)));hold on
plot(S_data);hold on
plot([1 204800],[Threshold Threshold]);
title('ȷ��ʧ��')


disp(['~���ٴ�ȷ��',char(fname),'�Ƿ���ʧ�ٹ�����'])
disp(['~��۲��Ƿ���Ƶ��ͼ�Ϲ۲쵽ʧ����תƵ�ʺ��䱶Ƶ�����ޣ���ʧ�ٹ���'])
disp(['~�Լ�ȷ��Thresholdֵ���Ƿ��ں���֮�䣿'])
prompt = '�۲�ʧ���˶����Ϊ1���۲�һ�ܱ仯Ϊ0��';
pause

[Pulse_RI,RI_Speed] = RI_RealTime(S_data,fs,Threshold); %ͨ�������źŻ�ȡʧ����λ��Ϣ 

% ÿȦ���Ҷ�Ӧ��ʧ����λ
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
[RI_order(:,1),RI_order(:,2)]=sort(+Phase);%һ�����ž�����ת����
    H_if_RI=RI_order(:,2)';
end     
%----------------%%%%---------------------%%%%%---------------------
%% 


solutime=1;
for H=H_if_RI  %ÿתһȦˢ��һ������һ��ʱ�䲽��%����RI_marker����
%% ����Ͱ������λ�ù�ϵ
distance_key_rotor=6;%����λ�ú����ҶƬ֮��ľ���
R_rotor=184.4002;%mm/��Ҷ�뾶
x=[2;4.41;7.68;11.14;14.39;17.41;20.73;23.77;26.63;34];
%y_2=[20.29;24.92;29.60;35.37;40.80;44.83;49.71;55.40;59.70;72.36];
y_2=[0;0;0;0;0;0;0;0;0;0];
distance_point_B1_R1=3.8*R_rotor*2*pi/29;
distance_point_B1_C1=5.3/360*R_rotor*2*pi;%17.0575mm;��C1�Ƶ�70mm����ֵΪ��70-17.0575��mm
distance_pointR1_R1=R_rotor*2*pi/29;
interval=Rotor_Speed(H)/60*R_rotor/fs*2*pi;
round_point=distance_pointR1_R1*29/interval;%һ�ܵĵ���
Len=floor(round_point)-2;%=130;%������ʾ�ڵ����;ȡһ��Χ��һ����ά��Ȧ
point_Pulse_rotor=round(distance_key_rotor/interval); %%%12000-��19*1.6092����13000-��12*1.7433����10000-��5*1.0740����11000-��-2*1.3598����125000-��16*1.6763��
point_rotor_rotor=round_point/29; %����ҶƬ֮��ľ������������һ��Ҷ���ĵ��� %144000/BPF
tran_C1=(70-distance_point_B1_C1)/interval;%C1ʵ�ʵ�λ�������λ��23.2235����λ��C1ǰ��23����λ,��24
point_B1_R1_interval=round(distance_point_B1_R1./interval);  %B1��R1��52.5761����λ��R1-R8ǰ��53����λ����54
point_B1_C1_interval=round((y_2(end)-y_2(1)-distance_point_B1_C1)./interval);  %%17.0575mm;��C1�Ƶ�70mm����ֵΪ��70-17.0575��mm
%% �����������
y_1=(1:Len+1)*interval;
X = repmat(x,size(y_1))';
yy_1 = repmat(y_1,size(x))';
yy_2=repmat(y_2,size(y_1))';
Y=yy_1+yy_2;
% 2ά���3ά��ת��
initial_y=Y(1,1);
Theta=(Y-initial_y)./R_rotor;
XX=R_rotor*cos(Theta);
YY=R_rotor*sin(Theta);
ZZ=X;

% �������ݣ�������ֵ
    rpm{solutime}(:,1) = signal(Pulse(H):Pulse(H)+Len,object(1));%������B1
    rpm{solutime}(:,2) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(2));
    rpm{solutime}(:,3) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(3));
    rpm{solutime}(:,4) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(4));
    rpm{solutime}(:,5) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(5));
    rpm{solutime}(:,6) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(6));
    rpm{solutime}(:,7) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(7));
    rpm{solutime}(:,8) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(8));
    rpm{solutime}(:,9) = signal((point_B1_R1_interval+Pulse(H)):(point_B1_R1_interval+Pulse(H)+Len),object(9));
    rpm{solutime}(:,10) = signal((point_B1_C1_interval+Pulse(H)):(point_B1_C1_interval+Pulse(H)+Len),object(10));%������C1
   % һȦ��29��ҶƬ��������л��֣���Ƕ�Ҷλ�ã�29�ݣ�:���
      for k=1:29  
        Pulse_solid(k)=round(1+point_rotor_rotor*(k-1));%�����Ϳ����ڱ�ת�ٲ���ҶƬ��
      end
     %rpm{solutime}(Pulse_solid+point_Pulse_rotor,:)=1.1*max(rpm{solutime}(:,end));  %�ڴ�Լpoint_rotor_rotor������Ѱ��ȷ��ҶƬλ�ã��ݶ��÷�����Ϊ����ָ��

 %%
     tsignal.surfaces(solutime).zonename='mysurface zone';
     tsignal.surfaces(solutime).x=XX;    %size 3x3 
     tsignal.surfaces(solutime).y=ZZ;    %size 3x3
     tsignal.surfaces(solutime).z=YY;    %size 3x3
     tsignal.surfaces(solutime).v(1,:,:)=rpm{solutime};%���ݼ����ź��ж�
     tsignal.surfaces(solutime).solutiontime=solutime;
     solutime=solutime+1;
end
 %% �������ݣ������ļ�
%title=''; 
rotorSpeed=floor(mean(Rotor_Speed)/100)*100;

if if_RI ==0
NAME = [strrep(char(fname),type,'-'),'-blade1-',date];  %�洢�ļ���
else
 NAME = ['DoubleLock-',strrep(char(fname),'.mat','-'),'-blade-',date];  %�洢�ļ���   
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
 
 
