function [Pulse,rotor_speed]=dianwoliu_keyRotation(thekeyData,fs)
   Tachometer = thekeyData(:,end);
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
    interval=floor((diff(the_Pulse))/29);
    for j=1:29
    for i=1:length(the_Pulse)-1
    L1_tip(i,j)=min(thekeyData((the_Pulse(i)+(j-2)*interval(i)):(the_Pulse(i)+(j-1)*interval(i)),1));
    L1_interval(i,j)=find(thekeyData((the_Pulse(i)+(j-2)*interval(i)):(the_Pulse(i)+(j-1)*interval(i)),1)==L1_tip(i,j))+((the_Pulse(i)+(j-2)*interval(i)))-1;
    L2_tip(i,j)=min(thekeyData((the_Pulse(i)+(j-2)*interval(i)):(the_Pulse(i)+(j-1)*interval(i)),2));
    L2_interval(i,j)=find(thekeyData((the_Pulse(i)+(j-2)*interval(i)):(the_Pulse(i)+(j-1)*interval(i)),2)==L2_tip(i,j))+((the_Pulse(i)+(j-2)*interval(i)))-1;
    
    end
    end
    figure
    plot(L1_tip(:,1:29));
    figure
    plot(diff(L1_interval(:,1:29)));
    figure
    plot(L2_tip(:,1:29));
    figure
    plot(diff(L2_interval(:,1:29)));
    figure
    L1_diff=diff(L1_interval');
    plot(L1_diff);
end
    
    
    
    
    