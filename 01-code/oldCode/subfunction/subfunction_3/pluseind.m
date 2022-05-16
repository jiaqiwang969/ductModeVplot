function I=pluseind(x)
%to caculate the pluse index of the signal
I=peakvalue(x)./meanamp(x);