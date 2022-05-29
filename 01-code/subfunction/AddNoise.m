function [Pnoise] = AddNoise(Pdirect,SNR,NumSM,NumMic,Nsnap)

% power_signal = sqrt(mean( (std(Pdirect',0,1) ).^2)); 
sigB = sqrt(mean(abs(Pdirect(:)).^2)/ (10^(0.1*SNR)));
% power_signal = mean( var(Pdirect',0,1) ); 
% power_noise = sqrt(power_signal/( 10^(0.1*SNR)));
Pnoise = Pdirect + sigB .*(randn(NumSM*NumMic,Nsnap)+1i*randn(NumSM*NumMic,Nsnap))/sqrt(2);
% NumSM:number of sequential measurements
% Y = randn(m,n) 或 Y = randn([m n])：返回一个m*n的随机项矩阵。
end