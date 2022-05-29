
function [tsignal2,Info]=computeDMD(rpm,tsignal,save_directory,name,rotor_Speed,fs,The_freq,Freq_dB,fntype)
for kk=1:length(rpm)
    len(kk)=length(rpm{1,kk});
end
[q1,q2]=min(len)
for kk=1:length(rpm)
    VORTALL(:,kk)=reshape(rpm{1,kk}(1:min(len),:),1,[])';%按列
end

X = VORTALL(:,1:end-1);
X2 = VORTALL(:,2:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r =30%length(rpm)-1;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi =real(X2*V*inv(S)*W);
EIGS=diag(eigs);
index_tecplot=find(imag(EIGS)>-1e-5); %tecplot只画上半平面
Info.rotor_Speed=rotor_Speed;
Info.dT=1/(rotor_Speed/60);
f=imag(log(EIGS)/Info.dT)/2/pi;

the_freq=[0:length(Phi) - 1]*204800/length(Phi)/(rotor_Speed/60);
data_freq=abs(fft(Phi(:,index_tecplot)))*2;

for k=1:29 
A_mode(:,k)=sum(data_freq(find(the_freq>k-0.5& the_freq<k+0.5),:));
end
Info.Phi=Phi(:,index_tecplot);
[c,v]=max(A_mode')
Info.EIGS=EIGS(index_tecplot);
Info.DMD_mode=[c' v' f(index_tecplot) v'*rotor_Speed/60-f(index_tecplot) v'*rotor_Speed/60+f(index_tecplot)];
Info.DMD_mode(:,6)=sign(sum(Freq_dB(round(Info.DMD_mode(:,4)/(The_freq(2)-The_freq(1)))+1,6:15)')- ... 
sum(Freq_dB(round(Info.DMD_mode(:,5)/(The_freq(2)-The_freq(1)))+1,6:15)'));
Info.DMD_mode(:,7)=(Info.DMD_mode(:,6)>-1).*Info.DMD_mode(:,4)+(Info.DMD_mode(:,6)==-1).*Info.DMD_mode(:,5);

tsignal2.Nvar= tsignal.Nvar;
tsignal2.varnames= tsignal.varnames;
for kk=1:length(index_tecplot)
tsignal2.surfaces(kk).zonename= tsignal.surfaces(q2(1)).zonename;
tsignal2.surfaces(kk).x= tsignal.surfaces(q2(1)).x;
tsignal2.surfaces(kk).y= tsignal.surfaces(q2(1)).y;
tsignal2.surfaces(kk).z= tsignal.surfaces(q2(1)).z;
tsignal2.surfaces(kk).v= reshape(Phi(:,index_tecplot(kk)),1,min(len),10);
%tsignal2.surfaces(kk).v(1,xuhao,:)=1.1*max(Phi(:,kk));
tsignal2.surfaces(kk).solutiontime=kk;
end
mat2tecplot(tsignal2,[save_directory,'\','TecplotDMD-IGV(',strrep(name,fntype,'-'),')'...
    num2str(round(rotor_Speed)),'rpm-',num2str(fs),'.plt']);
%% Plot DMD modes

% for i=10:2:20
%     figure
%     imagesc(reshape(real(Phi(:,i)),min(len),10)); % plot vorticity field
%     figure
%     imagesc(reshape(real(Phi(:,i)),min(len),10)); % plot vorticity field
% % 
% %     plotCylinder(reshape(real(Phi(:,i)),min(len),10),min(len),10);
% %     plotCylinder(reshape(imag(Phi(:,i)),min(len),10),min(len),10);
% end
% 
%%  Plot DMD spectrum
h=figure;set(gcf,'outerposition',get(0,'screensize'));%最大化
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
real_eigs=real(diag(eigs));
imag_eigs=imag(diag(eigs));
scatter(real_eigs,imag_eigs,'ok')
for i=1:length(index_tecplot)
text(real_eigs(index_tecplot(i)),imag_eigs(index_tecplot(i)),[num2str(i),'-',num2str(Info.DMD_mode(i,2))])
end
axis([-1.1 1.1 -1.1 1.1]);
saveas(h,[save_directory,'\','EigenValue-IGV(',strrep(name,fntype,'-'),')'...
    num2str(round(rotor_Speed)),'rpm-',num2str(fs),'.fig'])
saveas(h,[save_directory,'\','SpectrumAnalysis-IGV(',strrep(name,fntype,'-'),')'...
    num2str(round(rotor_Speed)),'rpm-',num2str(fs),'.png'])
% saveas(h,[save_directory,'\','Circle-',name,'.png'])
end