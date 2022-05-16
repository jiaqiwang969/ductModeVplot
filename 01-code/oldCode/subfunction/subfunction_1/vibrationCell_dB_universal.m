
function [Vib]=vibrationCell_dB_universal(signal,fs,object_1)

[Pulse,Rotor_Speed] = keyRotation_RealTime(signal(1:end,end),fs); %通过键向信号获取转速信息.去头去尾      
H_if_RI=1:(length(Rotor_Speed)-1);
solutime=1;
Len=min(diff(Pulse));
for H=H_if_RI  %每转一圈刷新一次网格，一个时间步！%按照RI_marker排序
for k=1:length(object_1)
    
    Vib{solutime}(:,k) = signal(Pulse(H):Pulse(H)+Len,object_1(k));%传感器B1~
end
    solutime=solutime+1;
end

end
 
 
