function [Pulse,Rotor_Speed]=keyRotation_RealTime(thekeyData,fs)
   Tachometer = thekeyData;
    Threshold = 4;
    Temp_Num =  find(Tachometer>Threshold); 
    k = 0;
    for i = 1:1:length(Temp_Num)-1
        if (Temp_Num(i+1) - Temp_Num(i)) > 1
            k = k+1;
            Pulse(k,1) = Temp_Num(i);
        end
    end
    
     Rotor_Speed= fs./((diff(Pulse)))*60;% rotor_speed每转一圈更新一次,并扩充成720000个点
     Rotor_Speed(end+1)=Rotor_Speed(end);
     Pulse(1)=[];Pulse(end)=[];   
     Rotor_Speed(1)=[];Rotor_Speed(end)=[];  %去头去尾,为了后面算相位
     %RI_marker=mod(Pulse,fs/f_RI)./(fs/f_RI)*360;    %假设Data序列一开始第一个数为RI的0相位
end
