function [Phase]=ModalWave_Phase(Pulse,Pulse_ModalWave)
% ÿȦ���Ҷ�Ӧ��maodl_wave��λ,����
for NumPhase=1:length(Pulse)
    temp_Phase=Pulse_ModalWave./Pulse(NumPhase);
    temp_order=find(temp_Phase<1);
     if isempty(temp_order)
        continue
     end
    if temp_order(end)==length(Pulse_ModalWave)
        break
    end
     Phase(NumPhase,:)=(Pulse(NumPhase)-Pulse_ModalWave(temp_order(end)))/(Pulse_ModalWave(temp_order(end)+1)-Pulse_ModalWave(temp_order(end)))*360;
end
end   