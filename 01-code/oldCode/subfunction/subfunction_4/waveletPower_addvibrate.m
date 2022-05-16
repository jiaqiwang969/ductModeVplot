%С���任����ʱ���źű任��ʱƵ��
%�ģ�������һ������fk_equalAngle��Ϊ����-20180409
function [PV_power]=waveletPower_addvibrate(sst1,fk_equalAngle)
%% waveletƵ������
variance = std(sst1).^2;
sst = (sst1 - mean(sst1))./sqrt(variance) ;

n = length(sst);
dt = 1/204800 ;
% round=[0:length(sst)-1]*dt*rotorSpeed_resample/60;  % construct round array
% xlim = [0,floor(round(end))];  % plotting range
pad = 1;      % pad the round series with zeroes (recommended)
dj = 0.125;    % this will do 4 sub-octaves per octave
s0 = 15/900*fk_equalAngle*dt;     % this says start at a scale of 6 months
j1 = 9/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.75;  % lag-1 autocorrelation for red noise background
mother = 'morlet';

% Wavelet transform:
for k=1:size(sst,2)
[wave,period,scale,coi] = wavelet(sst(:,k),dt,pad,dj,s0,j1,mother);
wave_band(k,:)=peakvalue(wave');%��ֵ��%��ȻҲ����������ʱ��ָ��
end
fk=1./period./(204800/fk_equalAngle);
PV_power = (abs(wave_band')).^2;        % compute wavelet power spectrum
end