function plot_A_error(mode_prop2,AA_re1,AA)
AA=AA.';
Y_real1=[real(AA),real(AA_re1)];
Y_imag1=[imag(AA),imag(AA_re1)];
Y_abs1=[abs(AA),abs(AA_re1)];
% Y_error_r1=abs(real(AA)-real(AA_re1))./real(AA);
% Y_error_i1=abs(imag(AA)-imag(AA_re1))./imag(AA);
Y_error_r1=real(AA)-real(AA_re1);
Y_error_i1=imag(AA)-imag(AA_re1);
% error_real=real(AA)-real(AA_re1);
% Y_error_r1=mean(error_real);
% error_imag=imag(AA)-imag(AA_re1);
% Y_error_i1=mean(error_imag);
Y_error_a1=abs(abs(AA)-abs(AA_re1))./abs(AA);
%%%%%%显示结果%%%%%%%%%
figure
subplot(4,1,1);
bar(Y_real1)
xlabel('mode order');ylabel('Mode coefficient (real part)');
legend('orignal','detected')
colormap(hot)

subplot(4,1,2);
bar(Y_imag1)
xlabel('mode order');ylabel('Mode coefficient(image part)');
legend('orignal','detected')
colormap(hot)

subplot(4,1,3);
bar(Y_error_r1)
xlabel('mode order');ylabel('error(real part)');
colormap(hot)

subplot(4,1,4);
bar(Y_error_i1)
xlabel('mode order');ylabel('error(image part)');
colormap(hot)

%%%%绘制3位柱状图%%%%
M_re1=[mode_prop2,AA_re1];
M=[mode_prop2,AA];
[y,M1_re1]=matrix_sort(M_re1);
figure
x=0:1:size(M1_re1,2);
x=x';
bar3(y,M1_re1)
xlabel('Radial mode n','Rotation',20);ylabel('Circumferential mode m','Rotation',-33);zlabel('Mode coefficient real part')
set(gca,'xticklabel',x)
set(gca,'ylim',[(y(1)-1) (y(end)+1)]);
title('detected mode')
[y,M1]=matrix_sort(M);
figure
bar3(y,M1)
xlabel('Radial mode n','Rotation',20);ylabel('Circumferential mode m','Rotation',-32);zlabel('Mode coefficient imag part')
set(gca,'xticklabel',x)
set(gca,'ylim',[(y(1)-1) (y(end)+1)]);
title('orignal mode')









