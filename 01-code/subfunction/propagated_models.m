%%%%扬声器（或传声器）位置
%   a——管道半径（m）
%   k—— 波数
function mode_prop2=propagated_models(k,Kappa)

% Kappa = Kappa(1:8,1:rings);%%%%%%首先根据圈数，挑选出最大径向模态阶数rings-1，
propa=sign(k^2-Kappa.^2);
[m,n]=find(propa==1);
mode_prop=[m-1,n-1];%%%根据截止频率条件，挑选出最大周向模态阶数
mode2=mode_prop(find(mode_prop(:,1)>0),:);%%挑选出m>0阶模态
mode2(:,1)=mode2(:,1)*-1;%%考虑-m阶模态
mode_prop2=[mode_prop;mode2];%%所有可能激发的模态数 m∈(-∞,+∞),n∈(0,+∞)