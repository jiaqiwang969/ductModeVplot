function [Pulse_RI,RI_Speed]=RI_RealTime(Tachometer,fs,Threshold) %R8

   
%    Tachometer =smooth(thekeyData,1000,'loess' );%��ͨ�˲�%��Ϊʧ�������ļ����ź�
%    
%    plot(Tachometer(1:10000)) 
   %Threshold = 3.08;%��ͷ�ģ���ѹֵ
    Temp_Num =  find(Tachometer<Threshold); 
    k = 0;
    for i = 1:1:length(Temp_Num)-1
        if (Temp_Num(i+1) - Temp_Num(i)) > 1
            k = k+1;
            Pulse_RI(k,1) = Temp_Num(i);
        end
    end
    
     RI_Speed= fs./((diff(Pulse_RI)))*60;% rotor_speedÿתһȦ����һ��,�������720000����
     RI_Speed(end+1)=RI_Speed(end);
     
     %RI_marker=mod(Pulse,fs/f_RI)./(fs/f_RI)*360;    %����Data����һ��ʼ��һ����ΪRI��0��λ
end
