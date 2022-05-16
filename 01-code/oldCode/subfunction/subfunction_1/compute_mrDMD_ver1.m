
function compute_mrDMD_ver1(rpm,tsignal,save_directory,name,average_degree,Pulse_blade,T4oneRound)
dt=T4oneRound/length(rpm);
T=T4oneRound;
for kk=1:length(rpm)
    len(kk)=length(rpm{1,kk});
end
[q1,q2]=min(len)
for kk=1:length(rpm)
    VORTALL(:,kk)=reshape(rpm{1,kk}(1:min(len),:),1,[])';%按列
end


%% compute mrDMD
L = 6; % number of levels
r = 10; % rank of truncation
mrdmd = mrDMD(VORTALL, dt, r, 10, L);
% compile visualization of multi-res mode amplitudes
[map, low_f] = mrDMD_map(mrdmd);
[L, J] = size(mrdmd);
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
imagesc(-map); 
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', floor(low_f*10)/10); 
set(gca, 'XTick', J/T*(0:T) + 0.5);
set(gca, 'XTickLabel', (get(gca, 'XTick')-0.5)/J*T);
axis xy;
xlabel('Time (sec)');
ylabel('Freq. (Hz)');
% colormap pink;
grid on;
colorbar('peer',axes1);

%% 直接在Matlab里作图
% for k=1:6
%     for kk=1:2^k/2
%         for kkk=1:size(mrdmd{k,kk}.Phi,2)
% h = figure;
% axes1 = axes('Parent',h);
% axis off
% hold(axes1,'on');
% surf(tsignal.surfaces(1).x,tsignal.surfaces(1).y,tsignal.surfaces(1).z,reshape(abs(mrdmd{1,1}.Phi(1:min(len)*10,1)),min(len),10),'Parent',axes1,'LineStyle','none');
% axis equal
% 
% title(['rmDMD-','[',num2str(k),',',num2str(kk),']-P',num2str(mrdmd{k,kk}.P(kkk)),'-rho',num2str(mrdmd{k,kk}.rho),'-omega',num2str(mrdmd{k,kk}.omega(kkk))])
% saveas(h,[save_directory,'\','rmDMD-','[',num2str(k),',',num2str(kk),']-P',num2str(mrdmd{k,kk}.P(kkk)),'-rho',num2str(mrdmd{k,kk}.rho),'-omega',num2str(mrdmd{k,kk}.omega(kkk)),'.png'])
% 
% pause(5*dt);
%         end
%     end
% end


%% Mesh_input
tsignal2.Nvar= tsignal.Nvar;
tsignal2.varnames= tsignal.varnames;
zhen=1;
for k=1:L
    for kk=1:2^k/2
       for kkk=1:size(mrdmd{k,kk}.Phi,2)%%

tsignal2.surfaces(zhen).zonename= tsignal.surfaces(q2(1)).zonename;
tsignal2.surfaces(zhen).x= tsignal.surfaces(q2(1)).x;
tsignal2.surfaces(zhen).y= tsignal.surfaces(q2(1)).y;
tsignal2.surfaces(zhen).z= tsignal.surfaces(q2(1)).z;
% tsignal2.surfaces(zhen).v= reshape(Xdmd(:,kk),1,min(len),10); 
% tsignal2.surfaces(zhen).v(1,Pulse_blade,:)=1.1*max(Phi(:,kk));
tsignal2.surfaces(zhen).v= reshape(abs(mrdmd{k,kk}.Phi(1:min(len)*10,kkk)),1,min(len),10);
%tsignal2.surfaces(zhen).v(1,Pulse_blade,:)=1.1*max(real(Xdmd(:,kk)));
tsignal2.surfaces(zhen).solutiontime=zhen;

zhen=zhen+1;
       end
    end
end

mat2tecplot(tsignal2,[save_directory,'\',name,'-mrDMD-d',num2str(average_degree),'r',num2str(r),'-reconstruction','.plt']);

% %%  Plot DMD spectrum
% h=figure
% theta = (0:1:100)*2*pi/100;
% plot(cos(theta),sin(theta),'k--') % plot unit circle
% hold on, grid on
% scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
% axis([-1.1 1.1 -1.1 1.1]);
% saveas(h,[save_directory,'\',name,'-d',num2str(average_degree),'r',num2str(r),'.png'])
end