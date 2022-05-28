%%%%构造系数矩阵G，用于模态识别%%%%
function [G,Kz]=matrix_G(mode_all,Kappa,k,a,mics_polar)
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
S=pi*a^2;
for i=1:row
    m_i=mode_all(i,1);n_i=mode_all(i,2);
    k_mn=Kappa(abs(mode_all(i,1))+1,abs(mode_all(i,2))+1);
    kz=sqrt(k^2-k_mn^2);%%%%轴向波数
    Kz(i)=sqrt(k^2-k_mn^2);
%     for j=1:size(mics_polar,1)
%         G(j,i)=besselj(abs(mode_all(i,1)),k_mn*mics_polar(j,1))*exp(-1j*mode_all(j,1)*mics_polar(j,2)*exp(1j*kz*mics_polar(j,3));
%     end
    G(:,count) = besselj(abs(mode_all(i,1)),k_mn*mics_polar(:,1)).*exp(-1j*mode_all(i,1).*mics_polar(:,2)).*exp(1j*kz*mics_polar(:,3));
    count=count+1;
end 






