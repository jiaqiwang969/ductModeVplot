
function [Vib]=vibrationCell_dB_universal(signal,fs,object_1)

[Pulse,Rotor_Speed] = keyRotation_RealTime(signal(1:end,end),fs); %ͨ�������źŻ�ȡת����Ϣ.ȥͷȥβ      
H_if_RI=1:(length(Rotor_Speed)-1);
solutime=1;
Len=min(diff(Pulse));
for H=H_if_RI  %ÿתһȦˢ��һ������һ��ʱ�䲽��%����RI_marker����
for k=1:length(object_1)
    
    Vib{solutime}(:,k) = signal(Pulse(H):Pulse(H)+Len,object_1(k));%������B1~
end
    solutime=solutime+1;
end

end
 
 
