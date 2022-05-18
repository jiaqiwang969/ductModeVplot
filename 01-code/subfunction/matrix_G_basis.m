%%%%构造系数矩阵G，用于模态识别%%%%
function [G,index_mn]=matrix_G_basis(f0,Radius,M,mic_loc,m,n)

z1 = (mic_loc(:,3).')/Radius;
theta1 = mic_loc(:,2).';


%% parameters
w=2*pi*f0/343*Radius;
rWall=1;

[Base] = BaseJ1(m,n);
beta=sqrt(1-M^2);
kappa_mn=sqrt(w^2-beta^2*Base.jmn_pm.^2);
kappa_sign=sign(w^2-beta^2*Base.jmn_pm.^2);

Eigp_mn=(-w*M-kappa_mn)/beta^2;  % right running


%Base.normValue(:,km).*
for km=1:length(m)
    Gmn1(:,:,km)=besselj(m(km),Base.jmn_pm(:,km)*rWall).*...
    exp(-i*(Eigp_mn(:,km)*z1)).*exp(i*m(km)*theta1);
end

G = reshape(permute(Gmn1,[2,1,3]),length(z1),n*length(m)); 
index_mn=[kron(m,ones(1,n));kron(ones(1,length(m)),[1:n])].';


% 删选cut-on模态，为了降低条件数
G(:,find(kappa_sign==-1))=[];
index_mn(find(kappa_sign==-1),:)=[];

