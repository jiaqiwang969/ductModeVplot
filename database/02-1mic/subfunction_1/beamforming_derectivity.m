function [rotor_speed]=beamforming_derectivity(signal,fs,Scale_data,R)
    the_freq = [0:length(signal)/2.56 - 1]*fs/length(signal);  %数据频域离散刻度
    temp_freq =fft(signal(:,:))*2./length(signal);
    temp_freq=temp_freq(1:length(signal)/2.56,:);
    [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
    Fk=[1/2,1,2,3];
     for K=1:4 
        initial=1;
     	for fai=linspace(-1/2*pi,1/2*pi,80)
        %转化为32个脉动球声源，然后再传播出去
        r1=0.05;%距离球声源的距离
        r_tube=0.2;%管道半径
        c=340;
        %通过相对坐标变换，确定测点相对声源的坐标
        x=R*cos(fai);y=R*sin(fai)-r_tube*cos(2*pi/32*1:32);z=r_tube*sin(2*pi/32*1:32);
        r2=sqrt((x).^2+(y).^2+(z).^2);
        A=r1*max(temp_freq(floor(rotor_speed/60*29*Fk(K)/(fs/length(signal)))+[floor(-2/(fs/length(signal))):floor(2/(fs/length(signal)))],1:32));
        p(initial)=abs(sum((A/r2).*exp(-i*2*pi*rotor_speed/60*29*Fk(K)/c*r2)));
        initial=initial+1;
    end
    subplot(2,2,K) 
    polar(linspace(-1/2*pi,1/2*pi,80),p,'r');hold on
    title([num2str(Fk(K)),'*BPF'],'FontSize',14);
     end
    suptitle({['-压气机进口指向性脉动球仿真'];['-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]})

end
   

   
  
 
