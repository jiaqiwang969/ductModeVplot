function Feature=features_wavelet(X,fk_equalAngle)
parfor i=1:length(X)  %≤¢––º∆À„
Feature{i}=waveletPower_addvibrate(X(i).v(1:fk_equalAngle,:),fk_equalAngle);
end