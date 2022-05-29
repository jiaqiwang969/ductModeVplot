function [Pulse,rotor_speed]=keyRotation(thekeyData,fs)
   Tachometer = thekeyData;
    Threshold = 4;
    Temp_Num =  find(Tachometer>Threshold); 
    k = 0;
    for i = 1:1:length(Temp_Num)-1
        if (Temp_Num(i+1) - Temp_Num(i)) > 1
            k = k+1;
            the_Pulse(k,1) = Temp_Num(i);
        end
    end
    
    fr = fs/(mean(diff(the_Pulse)));
    the_rotor_speed = fr*60;
    Pulse=the_Pulse;
    rotor_speed=floor(the_rotor_speed/10)*10;