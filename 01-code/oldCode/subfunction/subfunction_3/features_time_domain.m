function Feature=features_time_domain(X)

for i=1:length(X)
x=X(i).v;
%%ʱ������
x=x-mean(x);
PV=peakvalue(x)/10000; %��ֵ1
PPV=ppvalue(x)/10000;   %���ֵ2
AP=meanamp(x)/10000;   %ƽ����ֵ3
RP=rootamp(x)/10000;   %������ֵ4
RMS=rootmeansquare(x)/10000; %��Чֵ5
SK=skewness(x);   %���ָ��6
KU=kurtosis(x);  %�Ͷ�ָ��7
CF=peakind(x);   %��ֵָ��8
CLF=marginind(x); %ԣ��ָ��9
IF=pluseind(x);  %����ָ��10
SF=waveind(x);  %����ָ��11

Feature{i}=[PV;PPV;AP;RP;RMS;SK;KU;CF;CLF;IF;SF];
end
% figure 
% plot(A)



