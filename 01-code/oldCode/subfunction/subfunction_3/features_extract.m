function Feature=features(X)

parfor i=1:length(X)  %���м���
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
WE3=feature_WE3(x);%16��С����������12-28
SPEC=spectrumentropy(x);%��ֵ����28
ENV=Envelop(x);%��������29
Feature{i}=[PV;PPV;AP;RP;RMS;SK;KU;CF;CLF;IF;SF;WE3;SPEC;ENV];

end