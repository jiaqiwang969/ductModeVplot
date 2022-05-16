function [x,fs]=decimation_signal(old_x,old_fs,d_factor);

lnoldx=length(old_x);
lnx=floor(lnoldx/d_factor);
d=1:lnx;
l=d*d_factor;
x=old_x(l);
fs=old_fs/d_factor;

d=d/fs;
figure;
plot(d,x,'k');
showfft(x,fs);