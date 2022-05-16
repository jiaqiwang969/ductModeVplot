function [p_nm]=greenfunction(Boundary,m,Ratio,w,Tr,Omag,Mx,c02,rou0,k_nm,z_t,z_h,r0,x0)   
%numerical green's function calculation ,2018-10-10,wjq
for kk=1:length(k_nm)
xr = chebfun('xr', [Ratio, 1]);M_theta=Tr./xr+Omag*xr;
c02 = chebfun([fliplr(c02')]', [Ratio, 1]);rou0 = chebfun([fliplr(rou0')]', [Ratio,1]);
Mx = chebfun([fliplr(Mx')]',[Ratio,1]);
U_theta=2*M_theta/xr^2*diff(xr*M_theta)+M_theta^2/xr*(diff(rou0)/rou0-M_theta^2/xr/c02);
%% 计算格林函数
Phi=1-k_nm(kk)/w*Mx-m/w*M_theta/xr;Ri=(M_theta^2/xr)/c02*Phi+2*m/w*M_theta/xr^2;%(James-3.1.19)
A_rk=(U_theta-w^2*Phi^2)*Phi^2*w^2; %(James-3.5.21)
B_rk=w^2*Phi^2*((U_theta-w^2*Phi^2)*(1/xr-diff(rou0)/rou0)+diff(w^2*Phi^2-U_theta));%(James-3.5.22)
C_rkw2=(U_theta-w^2*Phi^2)^2*(Phi^2/c02-(k_nm(kk)/w)^2-(m/w)^2/xr^2)+Ri*(U_theta-w^2*Phi^2)*(Ri+Phi*(1/xr-diff(rou0)/rou0))+Ri*(diff(Phi*(w^2*Phi^2-U_theta)))-Phi*(w^2*Phi^2-U_theta)*(diff(Ri));%(James-3.1.21)
f1=(2*m/w*U_theta(1)/Phi(1)+U_theta(1)^2/c02(1))+(Boundary==2||Boundary==4)*sqrt(-1)*rou0(1)/z_t*(w*Phi(1)^2-U_theta(1)/w);%(James-3.1.22)&&(James-3.1.17)
f2=(2*m/w*U_theta(Ratio)/Phi(Ratio)/Ratio^2+U_theta(Ratio)^2/c02(Ratio))/Ratio-(Boundary==3||Boundary==4)*sqrt(-1)*rou0(1)/z_h*(w*Phi(Ratio)^2-U_theta(Ratio)/w);%(James-3.1.23)
L=chebop(Ratio, 1);
L.op=@(xr,g) A_rk*diff(g,2)+B_rk*diff(g,1)-w^2*C_rkw2.*g;%(James-3.5.20)
L1=L;L2=L;
L1.rbc=@(g) [g-1,diff(g)-f1];g1 =L1\0;
L2.lbc=@(g) [g-1,diff(g)-f2];g2 =L2\0;
% figure
% plot(real(g1),'.-');hold on;plot(real(g2),'.-');
%%
% df/dk
Phi_k=-Mx;Ri_k=-M_theta^2/xr/c02*Mx;
A_rk_k=(-w^2*2*Phi*Phi_k)*Phi^2*w^2+w^2*2*Phi*Phi_k*(U_theta-w^2*Phi^2);
B_rk_k=w^2*2*Phi*Phi_k*((U_theta-w^2*Phi^2)*(1/xr-diff(rou0)/rou0)+diff(w^2*Phi^2-U_theta))...
    +w^2*Phi^2*(-w^2*2*Phi*Phi_k*(1/xr-diff(rou0)/rou0)+diff(w^2*2*Phi*Phi_k-U_theta));
C_rk_kw2=2*(U_theta-w^2*Phi^2)*(-w^2)*2*Phi*Phi_k*(Phi^2/c02-(k_nm(kk)/w)^2-(m/w)^2/xr^2)...
    +(U_theta-w^2*Phi^2)^2*(2*Phi/c02*Phi_k-2*(k_nm(kk)/w))...
    +Ri_k*(U_theta-w^2*Phi^2)*(Ri+Phi*(1/xr-diff(rou0)/rou0))...
    +Ri*(-w^2*2*Phi*Phi_k)*(Ri+Phi*(1/xr-diff(rou0)/rou0))+Ri*(U_theta-w^2*Phi^2)*(Ri_k+(1/xr-diff(rou0)/rou0)*Phi_k)...
    +Ri_k*(diff(Phi*(w^2*Phi^2-U_theta)))+Ri*(diff(w^2*3*Phi^2*Phi_k-U_theta))...
    -(diff(Ri_k))*(Phi*(w^2*Phi^2-U_theta))-diff(Ri)*(w^2*3*Phi^2*Phi_k-U_theta);

b_k1=w^2*C_rk_kw2*g1-B_rk_k*(diff(g1))-A_rk_k*(diff(g1,2));
b_k2=w^2*C_rk_kw2*g2-B_rk_k*(diff(g2))-A_rk_k*(diff(g2,2));
f1_k=-2*m/w*U_theta(1)/Phi(1)^2*Phi_k(1)+(Boundary==2||Boundary==4)*sqrt(-1)*rou0(1)/z_t*w*2*Phi(1)*Phi_k(1);%(James-3.1.22)&&(James-3.1.17)
f2_k=-2*m/w*U_theta(Ratio)/Ratio^2/Phi(Ratio)^2*Phi_k(Ratio)-(Boundary==2||Boundary==4)*sqrt(-1)*rou0(Ratio)/z_h*w*2*Phi(Ratio)*Phi_k(Ratio);%(James-3.1.23)&&(James-3.1.17)
Lk1=L;Lk2=L;
Lk1.rbc=@(g) [g-0,diff(g)-f1_k];gk1 =Lk1\b_k1;
Lk2.lbc=@(g) [g-0,diff(g)-f2_k];gk2 =Lk2\b_k2;
% figure
% plot(real(gk1),'.-');hold on;plot(real(gk2),'.-');
W_k=gk1*diff(g2)+g1*diff(gk2)-g2*diff(gk1)-gk2*diff(g1); %(James-3.5.28)
kkk=1;
for xk0=x0
p_nm{kk,kkk}=sign(xk0)*((xr<r0)*sqrt(-1)*w/2/pi/r0/W_k(r0)*exp(sqrt(-1)*k_nm(kk)*xk0)*g1(r0)*g2+(xr>r0)*2*sqrt(-1)*w/4/pi/r0/W_k(r0)*exp(sqrt(-1)*k_nm(kk)*xk0)*g1*g2(r0));%(James-3.5.9)
kkk=kkk+1;
end
end