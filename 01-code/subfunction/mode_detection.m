%%%%本程序主要用于模态识别方法的仿真验证%%%%%%%%%%%%%
%%%%%模态识别，并分析误差
clc;clear all;
close all
addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\simulation\subprogram\'));
addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\'));
% addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\trial\subprogram\tikhonov\'));
% addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\trial\subprogram\TDMS\TDMS\'));
% addpath(genpath('E:\模态识别代码\mode_detect\mode_detect\trial\subprogram\TDMS\tdmsSubfunctions\'));
%% initial parameters 初始参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fa=4000; % 分析频率x
omega=2*pi*fa; 
c=340; k=omega/c; rho=1.225;pref=2e-5;
a=0.15; % 半径
S=pi*a^2; % 面积
rings=3;%%%%%%用了几圈扬声器
Rings=3;%用了几圈传声器
Ns=12;%单圈扬声器个数
load('Kappa.mat');
Kappa = Kappa/a; 
%% The propagation mode numbers 根据截至频率，计算出可传播模态%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode_prop2=propagated_models(k,Kappa); 
[row,col] = size(mode_prop2);
%% 扬声器阵列位置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('speaker.mat');%传声器安装位置，要求柱坐标
figure
xs = speaker(:,1).*cos(speaker(:,2));ys=speaker(:,1).*sin(speaker(:,2));
scatter3(xs,ys,speaker(:,3),'LineWidth',1.5);
%% 麦克风位置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('mic_jiaoda.mat');
mics=mic_jiaoda;
figure
xs = mics(:,1).*cos(mics(:,2));ys=mics(:,1).*sin(mics(:,2));
scatter3(xs,ys,mics(:,3),'LineWidth',1.5);
%% 令目标模态的模态系数为一个定值%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_p = 3; n_p = 0;%%%%%%需要产生的模态阶数
ind = ismember(mode_prop2,[m_p,n_p],'rows')';
posi = find(ind==1);
AA=zeros(1,row);AA(posi)=1;%%%%%目标模态系数为1，其他模态系数为0%%
%% 计算麦克风处的声压值%%%%%%%%%%%%%%%%%%%%%%%%%
[P1]=presure1(mics,mode_prop2,Kappa,k,a,AA);
P1=P1';
%% 计算传递矩阵G%%%%%%%%%%
[G,Kz]=matrix_G(mode_prop2,Kappa,k,a,mics);
%% 给添加噪声干扰%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR0=30;
P_real=noisegen(real(P1),SNR0);
P_imag=noisegen(imag(P1),SNR0);
P= P_real+1i*P_imag;
%% 方法:1：最小二乘法%%%%
cond(G)
AA_re1 = (G'*G)^-1*G'*P;%((G'*G)^-1)*G'*P2'
  %% 方法2：tikhonov正则化%%%%
% cond(G)
% [a,s,b] = csvd(G);
% lambda_l = l_curve (a,s,P1);
% AA_re1 = tikhonov (a,s,b,P1,lambda_l);
%% 方法:3：GMC%%%%
% lamda=0.01;
% zeta=10^(-0.3);
% rho=200*lamda*zeta;
% gamma=1.1;
% sc=1*10^-4;
% mlter=1000;
% [para_val,AA_re1,time]=ADMM_SRLASSO(G,P,lamda,rho,gamma,zeta, sc,mlter);
%% 方法:4：l2范数 %%%%
% lamda=0.1;
% cvx_begin
% variable x(36)
% minimize(0.5*(norm(G*x-P,'fro'))^2+lamda*norm(x,2))
% cvx_end
% AA_re1=x;
 %% 方法:4：l1范数 %%%%
% lamda=0.1;
% cvx_begin
% variable x(36)
% minimize(0.5*(norm(G*x-P,'fro'))^2+lamda*norm(x,1))
% cvx_end
% AA_re1=x;
%% 方法5：变分贝叶斯压缩感知
%      hyperpara.a0=1e-6;
%      hyperpara.b0=1e-6;
%      hyperpara.c0=1e-6;
%      hyperpara.d0=1e-6;
%      tol=1e-4;
%      plotflag=0;
%      Phi=G;
%      v=P;
%      niter=36;
%      [theta, posterior, bound]=BCSvb(Phi, v, hyperpara, niter, tol, plotflag)
%      AA_re1= theta;
 %% 方法6：贝叶斯压缩感知
 PHI=G;
 eta=1e-8;
 adaptive=0;
 optimal=1;
 scale=0.1;
 sigma2=SNR0;
 t=P;
q_re6= zeros(36,1);
[weights,used,sigma2,errbars] = BCS_fast_rvm(PHI,t,sigma2,eta,adaptive,optimal,scale)
 q_re6(used)= weights;
 %% 方法7：压缩感知  
% eta=1e-8;
% PHI=G;
% t=P;  
% a=SNR0;
% b=SNR0;
%  [weights,ML] = mt_CS(PHI,t,a,b,eta)
%  AA_re1=weights;
%% 计算误差与绘图
plot_A_error(mode_prop2,AA_re1,AA)
M_re1=[mode_prop2,AA_re1];
[y,M1_re1]=matrix_sort(M_re1);
% figure
% x=0:1:size(M1_re1,2);
% x=x';
% bar3(y,M1_re1)
% xlabel('Radial mode n','Rotation',20);ylabel('Circumferential mode m','Rotation',-33);zlabel('Mode coefficient real part')
% set(gca,'xticklabel',x)
% set(gca,'ylim',[(y(1)-1) (y(end)+1)]);
figure
abs_q1=20*log10(abs(AA_re1)/pref);%/(2*10-5)
abs_q1(find(abs_q1<0))=0;
bar(abs_q1)
xlabel('模态阶次');ylabel('模态系数A（dB）');
xlim([0 size(mode_prop2,1)+1])
colormap(hot)
ylim([max(abs_q1)-15 max(abs_q1)+5])
title('最小二乘法')
set(gca,'XTick',1:row);
set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(0,1)','(1,1)',...
    '(-1,0)','(-2,0)','(-3,0)','(-4,0)','(-1,1)'})
% AA1=AA'
% RMSE_real=zeros(1,60);
% RMSE_imag=zeros(1,60);
% RMSE_real(:,SNR0)=sqrt(mean((real(AA1(SNR0))-real(AA_re1)(SNR0)).^2));
% RMSE_imag(:,SNR0)=sqrt(mean((imag(AA1(SNR0))-imag(AA_re1)(SNR0)).^2));
figure
abs_q6=20*log10(abs(q_re6)/pref);%/(2*10-5)
abs_q6(find(abs_q6<0))=0;
bar(abs_q6)
xlabel('Modal order');ylabel('Modal coefficient A（dB）');
xlim([0 size(mode_prop2,1)+1])
colormap(hot)
ylim([max(abs_q6)-15 max(abs_q6)+5])
% title('贝叶斯压缩感知')
set(gca,'XTick',1:row);
set(gca,'XTickLabel',{'(0,0)','(1,0)','(2,0)','(3,0)','(4,0)','(5,0)','(6,0)','(7,0)','(8,0)',...
    '(0,1)','(1,1)','(2,1)','(3,1)','(4,1)','(0,2)','(1,2)','(-1,0)','(-2,0)','(-3,0)',...
    '(-4,0)','(-5,0)','(-6,0)','(-7,0)','(-8,0)','(-1,1)','(-2,1)','(-3,1)','(-4,1)','(-1,2)'})



 