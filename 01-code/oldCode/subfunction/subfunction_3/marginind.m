function L=marginind(x)
%to caculate the margin index of the signal
L=peakvalue(x)./rootamp(x);