%%%%此程序适用于：给定声源强度，求出管道内某位置处的声压值
%%%%%%2021.2.19%%%%%%%%%%%%%%%
function [P]=presure1(Location,mode_all,Kappa,k,a,AA)
%  Location——计算点处的空间位置，柱坐标系，m
%  mode_all--可产生的所有模态
%  omega——频率（rad/s）
%  rho——介质密度
%  k——波数
%  a——截面半径
%  AA——目标模态系数%%
%  P——声压计算值
N=size(Location,1);
row=size(mode_all,1);
S=pi*a^2;
for j=1:N
    for i=1:row
        m_i=mode_all(i,1);n_i=mode_all(i,2);
        k_mn=Kappa(abs(mode_all(i,1))+1,abs(mode_all(i,2))+1);
        kz=sqrt(k^2-k_mn^2);%%%%轴向波数
        pmn1(i)=AA(i)*besselj(abs(m_i),k_mn*Location(j,1))*exp(-1i*m_i*Location(j,2))*exp(1i*kz*Location(j,3));%%%%%%规定系数，直接求解声压
        clear ETA3
    end
    P(j)=sum(pmn1);%直接值
    clear pmn Amn pmn1
end














