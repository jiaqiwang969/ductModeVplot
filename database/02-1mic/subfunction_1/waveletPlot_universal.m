function [RotorSpeed]=waveletPlot_universal(signal,fs,Scale_data,object,objectName,testTime,fname,h)
    %function waveletPlot(RotorSpeed,location,fname,info_path,zfigure_wavelet_path)

fs_resample=20480;
[Pulse,RotorSpeed]=keyRotation(signal(:,end),fs);
sst=resample(signal(:,object),fs_resample,fs);
%variance = std(sst1)^2;
%sst = (sst1 - mean(sst1))/sqrt(variance) ;
n = length(sst);
dt = 1/fs_resample ;
round=[0:length(sst)-1]*dt*RotorSpeed/60;  % construct round array
xlim = [0,floor(round(end))];  % plotting range
pad = 1;      % pad the round series with zeroes (recommended)
dj = 0.125;    % this will do 4 sub-octaves per octave
s0 = 15*dt/(fs/fs_resample);     % this says start at a scale of 6 months
j1 = 9/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.75;  % lag-1 autocorrelation for red noise background
mother = 'morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum
fk=1./period/(RotorSpeed/60);
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(fk))):fix(log2(max(fk))));
contour(round,log2(fk),log2(power),log2(fk));  %*** or use 'contourfill'
%imagesc(round,log2(fs),log2(power));  %*** uncomment for 'image' plot
xlabel('Rotor revolutions (round)')
ylabel('Norm.Frequency（f/f_r_o_t）','FontSize',12)
title(' Wavelet Power Spectrum/lg2(Power)')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(fk),max(fk)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
grid on;



 end
   

% %waveletPlot(12000,location,char(fname(i_file)),info_path,zfigure_wavelet_path4paper)%$
% 
% %function waveletPlot(RotorSpeed,location,fname,info_path,zfigure_wavelet_path)
% fs_resample=20480;
% fs = 204800;   %￥采样频率Hz
% %load(fullfile(info_path,'kulite_transform_ab_0120.mat'));
% %Data = importdata(fullfile(location,fname));
% %Data=V2Pa(Data,info_path);
% %SST = Data(:,17);
% %figure
% %frequencyDomainPlote(SST,fs)
% %S_data =smooth(SST,1300,'loess' );%低通滤波%作为失速启动的键向信号
% %Threshold=mean(S_data);%默认值----->调整参数   
% %[Pulse,Rotor_Speed] = keyRotation_RealTime(Data(1:end,18),fs); %通过键向信号获取转速信息.去头去尾
% %[Pulse_ModalWave,ModalWave_Speed] = ModalWave_RealTime(S_data,fs,Threshold); %通过键向信号获取失速相位信息 
% %[Phase]=ModalWave_Phase(Pulse,Pulse_ModalWave);
% 
% %RotorSpeed=fix(mean(Rotor_Speed)/500)*500;
% %[the_freq,freq_dB]=frequencyDomainPlot_dB(SST,fs,10); %此处可放小波图作为对比
% sst1=resample(SST,fs_resample,fs);
% %S_data1=resample(S_data,fs_resample,fs);
% 
% %------------------------------------------------------ Computation
% 
% %variance = std(sst1)^2;
% %sst = (sst1 - mean(sst1))/sqrt(variance) ;
% 
% n = length(sst);
% dt = 1/fs_resample ;
% round=[0:length(sst)-1]*dt*RotorSpeed/60;  % construct round array
% xlim = [0,floor(round(end))];  % plotting range
% pad = 1;      % pad the round series with zeroes (recommended)
% dj = 0.125;    % this will do 4 sub-octaves per octave
% s0 = 15*dt/(fs/fs_resample);     % this says start at a scale of 6 months
% j1 = 9/dj;    % this says do 7 powers-of-two with dj sub-octaves each
% lag1 = 0.75;  % lag-1 autocorrelation for red noise background
% mother = 'morlet';
% 
% % Wavelet transform:
% [wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
% power = (abs(wave)).^2 ;        % compute wavelet power spectrum
% fk=1./period/(RotorSpeed/60);
% % Significance levels: (variance=1 for the normalized SST)
% [signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
% sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
% sig95 = power ./ sig95;         % where ratio > 1, power is significant
% 
% % Global wavelet spectrum & significance levels:
% for k_part=1:floor(round(end)/5)
% part{k_part}=find(round<5*k_part&round>5*(k_part-1));
% end
% for k_part=1:floor(round(end)/5)
% global_ws{k_part} =(sum(power(:,part{k_part})')./length(part{k_part}));   % round-average over all rounds
% end
% %global_ws{3} =(sum(power(:,part{3})')./length(part{3})); 
% dof = n - scale;  % the -scale corrects for padding at edges
% global_signif = wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);
% 
% % Scale-average between El Nino periods of 2--8 years
% avg = find((scale >= 2) & (scale < 8));
% Cdelta = 0.776;   % this is for the MORLET wavelet
% scale_avg = (scale')*(ones(1,n));  % expand scale --> (J+1)x(N) array
% scale_avg = power ./ scale_avg;   % [Eqn(24)]
% scale_avg = variance*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
% %scaleavg_signif = wave_signif(variance,dt,scale,2,lag1,-1,[0.002,0.0022],mother);
% 
% 
% %------------------------------------------------------ Plotting
% 
% h=figure('Visible', 'off');
% set(gcf,'outerposition',get(0,'screensize'));
% %--- Plot round series
% subplot('position',[0.1 0.75 0.65 0.2])
% plot(round,sst*sqrt(variance),'Color',[0 0.447058826684952 0.74117648601532]);hold on
% plot(round,S_data1);hold on
% plot([1 length(S_data1)],[Threshold Threshold]);
% set(gca,'Xgrid','on');set(gca,'gridlinestyle','--');
% set(gca,'XLim',xlim(:))
% xlabel('Rotor revolutions (round)')
% ylabel('Dynamic Pressure (Pa)')
% title(['a) Dynamic Pressure -',fname]);
% hold off
% 
% %--- Contour plot wavelet power spectrum
% fig=subplot('position',[0.1 0.11 0.65 0.53]);
% levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
% Yticks = 2.^(fix(log2(min(fk))):fix(log2(max(fk))));
% contour(round,log2(fk),log2(power),log2(fk));  %*** or use 'contourfill'
% %imagesc(round,log2(fs),log2(power));  %*** uncomment for 'image' plot
% xlabel('Rotor revolutions (round)')
% ylabel('Norm.Frequency（f/f_r_o_t）','FontSize',12)
% title('b) Wavelet Power Spectrum/lg2(Power)')
% set(gca,'XLim',xlim(:))
% set(gca,'YLim',log2([min(fk),max(fk)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks)
% grid on;
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% hold on
% contour(round,log2(period),sig95,[-99,1]);
% hold on
% % cone-of-influence, anything "below" is dubious
% plot(round,log2(coi))
% hold off
% % % 创建 textbox
% % annotation('textbox',...
% %     [0.675316455696202 0.462974310982475 0.0649789042513079 0.0392512083053589],...
% %     'String','Mode 1',...
% %      'HorizontalAlignment','center',...
% %     'FontSize',9,...
% %     'FitBoxToText','on',...
% %     'BackgroundColor',[1 1 1]);
% % 
% % % 创建 textbox
% % annotation('textbox',...
% %     [0.682469760900134 0.118167547697451 0.0590717311146893 0.0392512083053589],...
% %     'String','2*BPF',...
% %      'HorizontalAlignment','center',...
% %     'FontSize',9,...
% %     'FitBoxToText','on',...
% %     'BackgroundColor',[1 1 1]);
% % 
% % % 创建 textbox
% % annotation('textbox',...
% %     [0.682082981715886 0.179591521664968 0.0590717311146893 0.0392512083053589],...
% %     'String','1*BPF',...
% %      'HorizontalAlignment','center',...
% %     'FontSize',9,...
% %     'FitBoxToText','on',...
% %     'BackgroundColor',[1 1 1]);
% % 
% % % 创建 textbox
% % annotation('textbox',...
% %     [0.674874824191267 0.242595936865874 0.0666666680046275 0.0392512083053589],...
% %     'String','RI band',...
% %      'HorizontalAlignment','center',...
% %     'FontSize',9,...
% %     'FitBoxToText','on',...
% %     'BackgroundColor',[1 1 1]);
% % 
% % % 创建 textbox
% % annotation('textbox',...
% %     [0.65782137834036 0.521538438908625 0.0835443055378234 0.0392512083053589],...
% %     'String','1/2*Mode1',...
% %      'HorizontalAlignment','center',...
% %     'FontSize',9,...
% %     'FitBoxToText','on',...
% %     'BackgroundColor',[1 1 1]);
%     
% %--- Plot global wavelet spectrum
% subplot('position',[0.77 0.11 0.2 0.53])
% for k_part=[1:floor(round(end)/5)]
% plot(global_ws{k_part}/max(global_ws{k_part}),log2(fk),'DisplayName',[num2str((k_part-1)*5),'-',num2str((k_part)*5),'round'],'LineWidth',1); hold on
% end
% %plot(freq_dB/max(freq_dB),log2(the_freq),'DisplayName','dB','LineWidth',1); 
% hold off
% xlabel('Norm. Realtive Power (—)')
% title('c) Global Wavelet Spectrum')
% set(gca,'YLim',log2([min(fk),max(fk)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel','')
% set(gca,'XLim',[0,1.25])
% grid on;set(gca,'gridlinestyle','--');
% colorbar('peer',gca);
% legend1 = legend(gca,'off');
% set(legend1,...
%     'Position',[0.76926863456744 0.74909420289855 0.13115330636083 0.0612922705314016]);
% 
% saveas(h,[zfigure_wavelet_path,'\','compressor2stall-',fname,'WaveletPowerSpectrum-C2','.png']);
% close all
% end
