%%%%����ϵ������G������ģ̬ʶ��%%%%
function [G]=matrix_G_trial(mode_all,Kappa,k,a,mics_polar)
%  mode_all--�ɲ���������ģ̬
%  k��������
%  mics_polar����������λ��
%  omega����ƵƵ�ʣ�rad/s��
%  rho���������ܶ�
%  Kappa=Kappa(1:8,1:rings);%%%%%%���ȸ���Ȧ������ѡ�������ģ̬����rings-1��
%  a��������뾶
[row]=size(mode_all,1);
G=zeros(size(mics_polar,1),row);
count=1;
for i=1:row
    m_i=mode_all(i,1);n_i=mode_all(i,2);
    alpha_mn=Kappa(abs(mode_all(i,1))+1,abs(mode_all(i,2))+1);%alpha
    kz=sqrt(k^2-alpha_mn^2);%%%%������
    G(:,count) = besselj(abs(m_i),alpha_mn*mics_polar(:,1)).*exp(1j*m_i.*mics_polar(:,2)).*exp(-1j*kz*mics_polar(:,3));
    count=count+1;
end 






