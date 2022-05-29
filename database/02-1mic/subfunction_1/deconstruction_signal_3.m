function deconstruction_signal_3(signal,fs,object)
% fs=12000;
%data=resample(signal(:,object(1)),fs,fk);%qian
data=signal(:,object);%qian

% f=[4*10/fs 6*20/fs 2*350/fs 2*370/fs];
% A=[0 1 0];
% rp=0.153;
% rs=16.92;
% devp=1-10^(-rp/20);
% devs=10^(-rs/20);
% dev=[devp devs devp];
% [n,f0,A0,w]=remezord(f,A,dev);
% if rem(n,2)
%    n=n+1;
% end
% b=remez(n,f0,A0,w);
% freqz(b,1,length(b),1);
load filter_fs204800_f40_370.mat


data_filter = filter(b,1,data);


N = length(data_filter);
x = (data_filter - mean(data_filter))/std(data_filter);
xe= (data - mean(data))/std(data);

[row,col] = size(x);
if row<col
    x = x';
end


% figure
% subplot(212)
% [psd_x,f] = pwelch(x,[],[],round(N/4),fs);
% 
% plot(f,psd_x)
% hold on
% plot(f1,psd_x1)


% Demodulation
data_fault_filter=x;
Nfft = length(data_fault_filter);
f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
hilbert_data_fault_filter = hilbert(data_fault_filter);
envelope_spec_filter = abs(fft(abs(hilbert_data_fault_filter)));
% figure
createfigure_paper2_envelopesingal(f,[envelope_spec_filter(1:length(f))./Nfft]')
title('ÂË²¨Ö®ºó')
%plot(f,envelope_spec_filter(1:Nfft/2)./Nfft)
% hold on
% e1=xe;
% Nfft = length(e1);
% f = [0:fs/Nfft:(Nfft-1)*fs/(2*Nfft)];
% hilbert_e1 = hilbert(e1);
% envelope_spec_e1 = abs(fft(abs(hilbert_e1)));
% %plot(f,envelope_spec_e1(1:Nfft/2)./Nfft,'r')
% createfigure_paper2_envelopesingal(f,envelope_spec_e1(1:length(f))./Nfft)
% title('no lvbo')
% 
% figure 
% plot(e1)
% hold on 
% plot(abs(hilbert_e1))
% title('no lvbo')


% figure 
% plot(data_fault_filter)
% hold on 
% plot(abs(hilbert_data_fault_filter))
% title('lvbo')


end
 
      


   
