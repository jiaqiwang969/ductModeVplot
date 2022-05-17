%%%%构造系数矩阵G，用于模态识别%%%%
function [G]=matrix_G_trial(mode_all,Kappa,k,a,mics_polar)
%  mode_all--可产生的所有模态
%  k——波数
%  mics_polar——传声器位置
%  omega——频频率（rad/s）
%  rho——介质密度
%  Kappa=Kappa(1:8,1:rings);%%%%%%首先根据圈数，挑选出最大径向模态阶数rings-1，
%  a——截面半径
[row]=size(mode_all,1);
G=zeros(size(mics_polar,1),row);
count=1;
for i=1:row
    m_i=mode_all(i,1);n_i=mode_all(i,2);
    alpha_mn=Kappa(abs(mode_all(i,1))+1,abs(mode_all(i,2))+1);%alpha
    kz=sqrt(k^2-alpha_mn^2);%%%%轴向波数
    G(:,count) = besselj(abs(m_i),alpha_mn*mics_polar(:,1)).*exp(1j*m_i.*mics_polar(:,2)).*exp(-1j*kz*mics_polar(:,3));
    count=count+1;
end 






