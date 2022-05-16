%WAVETEST Example Matlab script for WAVELET, using NINO3 SST dataset
%
% See "http://paos.colorado.edu/research/wavelets/"
% Written January 1998 by C. Torrence
%
% Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
%   changed all "log" to "log2", changed logarithmic axis on GWS to
%   a normal axis.
%close all
clc
clear
fs=204800;
fk=20480;
RotorSpeed=10000;
info_path='H:\compressor2stall-2018-1-20\动画预测一体化程序\code_info';
load(fullfile(info_path,'kulite_transform_ab_0120.mat'));
[fname,location]=uigetfile({'*.*';'*.mat'},'r');%选择文件
testdata = importdata(fullfile(location,fname));
SST = testdata(:,4)*kulite_transform_ab(1,1)+kulite_transform_ab(1,2);
sst=resample(SST,fk,fs);
sst1=sst;

%------------------------------------------------------ Computation

% normalize by standard deviation (not necessary, but makes it easier
% to compare with plot on Interactive Wavelet page, at
% "http://paos.colorado.edu/research/wavelets/plot/"
variance = std(sst1)^2;
sst = (sst1 - mean(sst1))/sqrt(variance) ;

n = length(sst);
dt = 1/fk ;
round=[0:length(sst)-1]*dt*RotorSpeed/60;  % construct round array
xlim = [0,floor(round(end))];  % plotting range
pad = 1;      % pad the round series with zeroes (recommended)
dj = 0.125;    % this will do 4 sub-octaves per octave
s0 = 15*dt/(fs/fk);    % this says start at a scale of 6 months
j1 = 9/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.75;  % lag-1 autocorrelation for red noise background
mother = 'morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum
fk=1./period/(RotorSpeed/60);
% Significance levels: (variance=1 for the normalized SST)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:
part{1}=find(round<15);part{2}=find(round>20&round<35);%part{3}=find(round>20&round<30);

global_ws{1} =(sum(power(:,part{1})')./length(part{1}));   % round-average over all rounds
global_ws{2} =(sum(power(:,part{2})')./length(part{2})); 
%global_ws{3} =(sum(power(:,part{3})')./length(part{3})); 
dof = n - scale;  % the -scale corrects for padding at edges
global_signif = wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);

% Scale-average between El Nino periods of 2--8 years
avg = find((scale >= 2) & (scale < 8));
Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scale')*(ones(1,n));  % expand scale --> (J+1)x(N) array
scale_avg = power ./ scale_avg;   % [Eqn(24)]
scale_avg = variance*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
%scaleavg_signif = wave_signif(variance,dt,scale,2,lag1,-1,[0.002,0.0022],mother);

whos

%------------------------------------------------------ Plotting

%--- Plot round series
subplot('position',[0.1 0.75 0.65 0.2])
plot(round,sst1-mean(sst1))
set(gca,'Xgrid','on');set(gca,'gridlinestyle','--');
set(gca,'XLim',xlim(:))
xlabel('Rotor revolutions (round)')
ylabel('Dynamic Pressure (Pa)')
title(['a) Compressor2Stall Dynamic Pressure (',num2str(RotorSpeed),'rpm)']);
hold off

%--- Contour plot wavelet power spectrum
fig=subplot('position',[0.1 0.11 0.65 0.53])
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(fk))):fix(log2(max(fk))));
contour(round,log2(fk),log2(power),log2(fk));  %*** or use 'contourfill'
%imagesc(round,log2(fs),log2(power));  %*** uncomment for 'image' plot
xlabel('Rotor revolutions (round)')
ylabel('Norm.Frequency（f/f_r_o_t）','FontSize',12)
title('b) Wavelet Power Spectrum/lg2(Power)')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(fk),max(fk)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
grid on;
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(round,log2(period),sig95,[-99,1]);
hold on
% cone-of-influence, anything "below" is dubious
plot(round,log2(coi))
hold off
% 创建 textbox
annotation('textbox',...
    [0.675316455696202 0.462974310982475 0.0649789042513079 0.0392512083053589],...
    'String','Mode 1',...
     'HorizontalAlignment','center',...
    'FontSize',9,...
    'FitBoxToText','on',...
    'BackgroundColor',[1 1 1]);

% 创建 textbox
annotation('textbox',...
    [0.682469760900134 0.118167547697451 0.0590717311146893 0.0392512083053589],...
    'String','2*BPF',...
     'HorizontalAlignment','center',...
    'FontSize',9,...
    'FitBoxToText','on',...
    'BackgroundColor',[1 1 1]);

% 创建 textbox
annotation('textbox',...
    [0.682082981715886 0.179591521664968 0.0590717311146893 0.0392512083053589],...
    'String','1*BPF',...
     'HorizontalAlignment','center',...
    'FontSize',9,...
    'FitBoxToText','on',...
    'BackgroundColor',[1 1 1]);

% 创建 textbox
annotation('textbox',...
    [0.674874824191267 0.242595936865874 0.0666666680046275 0.0392512083053589],...
    'String','RI band',...
     'HorizontalAlignment','center',...
    'FontSize',9,...
    'FitBoxToText','on',...
    'BackgroundColor',[1 1 1]);

% 创建 textbox
annotation('textbox',...
    [0.65782137834036 0.521538438908625 0.0835443055378234 0.0392512083053589],...
    'String','1/2*Mode1',...
     'HorizontalAlignment','center',...
    'FontSize',9,...
    'FitBoxToText','on',...
    'BackgroundColor',[1 1 1]);
    
%--- Plot global wavelet spectrum
subplot('position',[0.77 0.11 0.2 0.53])

plot(global_ws{1}/max(global_ws{1}),log2(fk),'DisplayName','0-15 round','LineWidth',3); hold on
plot(global_ws{2}/max(global_ws{2}),log2(fk),'DisplayName','15-35 round','LineWidth',3,'Color',[1 0.843137264251709 0]); hold on

hold off
xlabel('Norm. Realtive Power (—)')
title('c) Global Wavelet Spectrum')
set(gca,'YLim',log2([min(fk),max(fk)]), ...
	'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel','')
set(gca,'XLim',[0,1.25])
grid on;set(gca,'gridlinestyle','--');
colorbar('peer',gca);
legend1 = legend(gca,'show');
set(legend1,...
    'Position',[0.76926863456744 0.74909420289855 0.13115330636083 0.0612922705314016]);



