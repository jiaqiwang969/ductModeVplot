function [Pulse_RI,Pulse_RI_low1,Pulse_RI_high1,RI_Speed]=RI_PhaseLockPoint(Tachometer,fs,Threshold) %R8

   
%    Tachometer =smooth(thekeyData,1000,'loess' );%低通滤波%作为失速启动的键向信号
%    
%    plot(Tachometer(1:10000)) 
   %Threshold = 3.08;%回头改，电压值
    Temp_Num =  find(Tachometer<Threshold); 
    k = 0;
    for i = 1:1:length(Temp_Num)-1
        if (Temp_Num(i+1) - Temp_Num(i)) > 1
            k = k+1;
            Pulse_RI(k,1) = Temp_Num(i);
        end
    end
    
  
    RI_Speed= fs./((diff(Pulse_RI)))*60;% rotor_speed每转一圈更新一次,并扩充成720000个点
    RI_Speed(end+1)=RI_Speed(end);
    %RI_marker=mod(Pulse,fs/f_RI)./(fs/f_RI)*360;    %假设Data序列一开始第一个数为RI的0相位
   
    for kk=1:length(Pulse_RI)-1
        [Pulse_RI_low2(kk,1),temp1]=min(Tachometer(Pulse_RI(kk):Pulse_RI(kk+1)));
        Pulse_RI_low1(kk,1)=temp1+Pulse_RI(kk); 
        [Pulse_RI_high2(kk,1),temp2]=max(Tachometer(Pulse_RI(kk):Pulse_RI(kk+1)));
        Pulse_RI_high1(kk,1)=temp2+Pulse_RI(kk); 
    end    
