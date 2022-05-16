clear;
clc;
close all;
set(0,'DefaultAxesFontSize',14,'DefaultAxesLineWidth',1,'DefaultFigurePosition',[100 100 780 360]);
fs=25600;%����Ƶ��
fr=12;%תƵ
%fc=12;%���ּ�
fo=37;%��Ȧͨ��Ƶ��
fb=143;%���������Ƶ��
fi=54;%��Ȧͨ��Ƶ
%���ݳ���
%fp=fi;%������Ȧ
%Q=fr;
fp=fo;%������Ȧ
Q=0;
fn=2000;
N=40;
A=1;
var_tao=0.00;
x=create_bearing_signal(N,A,fs,fn,fp,var_tao,Q);
Num=length(x);
df=fs/Num;
figure;
plot(x);
Nt=1024;
Ndt=32;
Ntau=256;
alpha=fo;
cyc_crocorr_certainalpha_gbi(x',alpha,fs,Nt,Ndt,Ntau);
% for i=1:floor(Num/5.12)
%     alpha=i*df;
%     R=cyc_crocorr_certainalpha_gbi(x',alpha,fs,Nt,Ndt,Ntau);
%     d(i)=degree_of_cyclostationarity(R,fs);
% end
% f=(1:floor(Num/5.12))*df;
% figure;
% plot(f,d);