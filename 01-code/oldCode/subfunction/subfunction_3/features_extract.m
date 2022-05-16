function Feature=features(X)

parfor i=1:length(X)  %并行计算
x=X(i).v;
%%时域特征
x=x-mean(x);
PV=peakvalue(x)/10000; %峰值1
PPV=ppvalue(x)/10000;   %峰峰值2
AP=meanamp(x)/10000;   %平均幅值3
RP=rootamp(x)/10000;   %方根幅值4
RMS=rootmeansquare(x)/10000; %有效值5
SK=skewness(x);   %歪度指标6
KU=kurtosis(x);  %峭度指标7
CF=peakind(x);   %峰值指标8
CLF=marginind(x); %裕度指标9
IF=pluseind(x);  %脉冲指标10
SF=waveind(x);  %波形指标11
WE3=feature_WE3(x);%16个小波包能量熵12-28
SPEC=spectrumentropy(x);%幅值谱熵28
ENV=Envelop(x);%包络谱熵29
Feature{i}=[PV;PPV;AP;RP;RMS;SK;KU;CF;CLF;IF;SF;WE3;SPEC;ENV];

end