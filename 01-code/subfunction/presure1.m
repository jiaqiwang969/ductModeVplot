%%%%�˳��������ڣ�������Դǿ�ȣ�����ܵ���ĳλ�ô�����ѹֵ
%%%%%%2021.2.19%%%%%%%%%%%%%%%
function [P]=presure1(Location,mode_all,Kappa,k,a,AA)
%  Location��������㴦�Ŀռ�λ�ã�������ϵ��m
%  mode_all--�ɲ���������ģ̬
%  omega����Ƶ�ʣ�rad/s��
%  rho���������ܶ�
%  k��������
%  a��������뾶
%  AA����Ŀ��ģ̬ϵ��%%
%  P������ѹ����ֵ
N=size(Location,1);
row=size(mode_all,1);
S=pi*a^2;
for j=1:N
    for i=1:row
        m_i=mode_all(i,1);n_i=mode_all(i,2);
        k_mn=Kappa(abs(mode_all(i,1))+1,abs(mode_all(i,2))+1);
        kz=sqrt(k^2-k_mn^2);%%%%������
        pmn1(i)=AA(i)*besselj(abs(m_i),k_mn*Location(j,1))*exp(-1i*m_i*Location(j,2))*exp(1i*kz*Location(j,3));%%%%%%�涨ϵ����ֱ�������ѹ
        clear ETA3
    end
    P(j)=sum(pmn1);%ֱ��ֵ
    clear pmn Amn pmn1
end














