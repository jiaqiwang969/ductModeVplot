function [RPM1,tsignal,Pulse_blade,T,signal_X]=pltPlot_dB_universal_ver1(signal,fs,object,objectName,testTime,fname,save_directory,if_RI,type,average_degree)

[Pulse,Rotor_Speed] = keyRotation_RealTime(signal(1:end,end),fs); %ͨ�������źŻ�ȡת����Ϣ.ȥͷȥβ      
N_blade=4;%���һ������-3��4��Ҷ��Ϊһ��-������������
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

%% ���ʧ��˫���๦��
%---------------%%-%%%%--------��ע�͵�-------------%%%%%---------------------
if if_RI ==1

figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes;
hold(axes1,'on');
subplot(4,1,1)
S_data =smooth(signal(:,4),200,'loess' );%��ͨ�˲�%��Ϊʧ�������ļ����ź�
Threshold=mean(S_data);%Ĭ��ֵ----->�������� 
plot1=plot(signal(:,4));
set(plot1(1),'Color',[0.862745106220245 0.862745106220245 0.862745106220245]);
hold on
plot(S_data,'LineWidth',2);hold on
plot([1 204800],[Threshold Threshold],'LineStyle','--');


axes1 = axes;
hold(axes1,'on');
subplot(4,1,2)
S_data =smooth(signal(:,5),200,'loess' );%��ͨ�˲�%��Ϊʧ�������ļ����ź�
Threshold=mean(S_data);%Ĭ��ֵ----->�������� 
plot1=plot(signal(:,5));
set(plot1(1),'Color',[0.862745106220245 0.862745106220245 0.862745106220245]);
hold on
plot(S_data,'LineWidth',2);hold on
plot([1 204800],[Threshold Threshold],'LineStyle','--');

axes1 = axes;
hold(axes1,'on');
subplot(4,1,3)
S_data =smooth(signal(:,6),200,'loess' );%��ͨ�˲�%��Ϊʧ�������ļ����ź�
Threshold=mean(S_data);%Ĭ��ֵ----->�������� 
plot1=plot(signal(:,6));
set(plot1(1),'Color',[0.862745106220245 0.862745106220245 0.862745106220245]);
hold on
plot(S_data,'LineWidth',2);hold on
plot([1 204800],[Threshold Threshold],'LineStyle','--');


axes1 = axes;
hold(axes1,'on');
subplot(4,1,4)
S_data =smooth(signal(:,7),200,'loess' );%��ͨ�˲�%��Ϊʧ�������ļ����ź�
Threshold=mean(S_data);%Ĭ��ֵ----->�������� 
plot1=plot(signal(:,7));
set(plot1(1),'Color',[0.862745106220245 0.862745106220245 0.862745106220245]);
hold on
plot(S_data,'LineWidth',2);hold on
plot([1 204800],[Threshold Threshold],'LineStyle','--');



end
 
      


   
