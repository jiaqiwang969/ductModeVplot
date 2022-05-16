function pltPlot_dB_universal_ver1(signal,fs,object,objectName,testTime,fname,save_directory,if_RI,type,average_degree)
% secondPaper-双锁相低通滤波参数讨论discuss1-
%---》S_data =smooth(signal(:,object(9)),900,'loess' );%低通滤波%作为失速启动的键向信号
%---》

    
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
    span=500:100:5000;
    for k_object=1:length(object) 
        for k_span=1:length(span) 
        S_data{k_object,k_span} =smooth(signal(:,object(k_object)),span(k_span),'loess' );%低通滤波%作为失速启动的键向信号
        Threshold{k_object,k_span}=mean(S_data{k_object,k_span});%默认值----->调整参数 
        [Pulse_RI{k_object,k_span},Pulse_RI_low1{k_object,k_span},Pulse_RI_high1{k_object,k_span},RI_Speed] = RI_PhaseLockPoint(S_data{k_object,k_span},fs,Threshold{k_object,k_span}); %锁相点
        evenness{k_object}(k_span,1)=max(diff(Pulse_RI{k_object,k_span}))-min(diff(Pulse_RI{k_object,k_span}));
        evenness{k_object}(k_span,2)=max(diff(Pulse_RI_low1{k_object,k_span}))-min(diff(Pulse_RI_low1{k_object,k_span}));
        evenness{k_object}(k_span,3)=max(diff(Pulse_RI_high1{k_object,k_span}))-min(diff(Pulse_RI_high1{k_object,k_span}));
  
        end
       evenness_min=min(min(evenness{k_object}));
       [span_opt,column]=find(evenness{k_object}==evenness_min);
       createfigure(signal,object,span,k_object,span_opt,S_data,Threshold,column,evenness_min,Pulse_RI,Pulse_RI_low1,Pulse_RI_high1,save_directory);      
    end            
end
end 
function createfigure(signal,object,span,k_object,span_opt,S_data,Threshold,column,evenness_min,Pulse_RI,Pulse_RI_low1,Pulse_RI_high1,save_directory) 
opt_LockMethod=['zero';'lowe';'high'];
objectName1=['B-1';'R-1';'R-2';'R-3';'R-4';'R-5';'R-6';'R-7';'R-8';'C-1';];
h=figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes;
hold(axes1,'on');
plot1=plot(signal(:,object(k_object)));
set(plot1(1),'Color',[0.862745106220245 0.862745106220245 0.862745106220245]);
hold on
plot(S_data{k_object,span_opt},'LineWidth',3);hold on
plot([1 204800],[Threshold{k_object,span_opt} Threshold{k_object,span_opt}],'LineStyle','--');
title([objectName1(k_object,:),'optLockMethod->',opt_LockMethod(column,:),'-',num2str(span(span_opt)),'-',num2str(evenness_min)]);
hold on
plot(Pulse_RI{k_object,span_opt},S_data{k_object,span_opt}(Pulse_RI{k_object,span_opt}),'MarkerSize',25,'Marker','.','LineWidth',15,'LineStyle','none')
hold on
plot(Pulse_RI_low1{k_object,span_opt},S_data{k_object,span_opt}(Pulse_RI_low1{k_object,span_opt}),'MarkerSize',25,'Marker','.','LineStyle','none','Color',[1 0 0])
hold on
plot(Pulse_RI_high1{k_object,span_opt},S_data{k_object,span_opt}(Pulse_RI_high1{k_object,span_opt}),'MarkerSize',25,'Marker','.','LineStyle','none','Color',[1 1 0])

xlim(axes1,[Pulse_RI{k_object,span_opt}(2) Pulse_RI{k_object,span_opt}(6)]);
saveas(h,[save_directory,'\','PLA-discuss1-',objectName1(k_object,:),'optLockMethod-',opt_LockMethod(column,:),'-',num2str(span(span_opt)),'-',num2str(evenness_min),'.png'])
saveas(h,[save_directory,'\','PLA-discuss1-',objectName1(k_object,:),'optLockMethod-',opt_LockMethod(column,:),'-',num2str(span(span_opt)),'-',num2str(evenness_min),'.fig'])

end
   
