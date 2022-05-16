%%����������ģ̬��Ƶ�ʽ�ֹ��ϵ��ͼ����������2019-9-26�޸�;
%������FFT�����ļ���Աȣ�����μ�wavemode_calculation_FFT.m����������2019-09-30
%wjq, 2019-10-6�޸ģ��������й���
%��ת��ϻģ̬����
clc
clear
% close all
i_cpsd_method = 1;  %CPSD���������1:ʱ���ף�2���ױ���
i_spectrum_method = 2;    %1:fftδ�ֶ�ƽ����2��fft�ֶ�ƽ����


subfunction_path1='H:\����Ԥ��һ�廯����ͨ�ð�_��ת��ϻ_ver1\subfunction\subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.txt';'*.*'},'mat�����ļ���ȡ','MultiSelect','on');%MultiSelect��ѡ
load([location,'\','����˵��','\','parameter12000.mat']); %ѡ���ļ���������
disp(Note);
% % //======����ͼ����ָ���ļ���===============  
save_directory = [strrep(location,'Database','11-16-����������ģ̬(CPSD)'),date];  %Ƶ��ͼ�洢�ļ���
   
    
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('�ļ��д��ڣ�');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
for i_file=1:length(fname)
    close all
    Data = importdata(fullfile(location,char(fname(i_file)))); %ѡ���ļ���������
    Data=V2Pa_Universal(Data,kulite_transform_ab);
    Data(:,1:end-1)=Data(:,1:end-1)-mean(Data(:,1:end-1));
    DATA{i_file}=Data(:,[1:13 14]);
% %%  1-Ƶ��ͼ    
%     [freq,rotor_speed,freq_dB]=spectrum_calculation_FFT(Data,fs,17,i_spectrum_method,save_directory,fname,i_file);
end
%    
% %% 2-ģ̬ͼ-by CPSD      
%    
% [SPECTRA,the_freq]=wavemode_calculation_CPSD_combine(DATA,fs,12,save_directory,fname,16000);
  L_signal = fs*5;%1������
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft = 2^(ceil(log2(L_seg))+1);  
    nk=12
    rotor_speed=6000;
     for i_file=1:length(fname) 
%           phase_info{i_file}=cpsd(DATA{i_file}(:,13),DATA{1}(:,13),Wind,Noverlap,Nfft,fs)./...
%                 abs(cpsd(DATA{i_file}(:,13),DATA{1}(:,13),Wind,Noverlap,Nfft,fs)); 
          for k=1:nk
            [temp,freq] = cpsd(DATA{i_file}(:,k),DATA{i_file}(:,13),Wind,Noverlap,Nfft,fs);
             CC1(:,i_file+length(fname)*(k-1)) = temp;%./phase_info{i_file};
        end
     end
    df=freq(2)-freq(1);
    
    nk_enlarge=length(fname)*nk;
    m=-nk_enlarge/2:nk_enlarge/2;
    for k=1:length(m)
        a_mf(:,k)=1/nk_enlarge*CC1(:,1:nk_enlarge)*exp(m(k)*i*2*pi*(1:nk_enlarge)/nk_enlarge).'; 
    end
    
      h=figure('Visible', 'on'); set(gcf,'outerposition',get(0,'screensize'));%���
     GAMMA =10*log10(abs(a_mf)/4e-10);
%      GAMMA(find(GAMMA<75))=50;
     imagesc([-nk_enlarge/2:nk_enlarge/2],freq,GAMMA);ylim([1,rotor_speed/60*29*3.2]); xlim([-20,20]); axis xy;
     title(['��ֹģ̬����',' -CPSD method ', fname{1, 1}])


% f = 3000;
% %% ---------��Դ����------------------------------------------------------------
% row = 1.29;                               %�ܶ�
speed = 344;                              %����
% a = 0.01;                                 %��Դ�뾶
% v0 = 2.5e-4;                              %��Դ�����ٶ�
zR = 0.05;                                 %�ؽ�����
zH = 3.05;                                 %��������
% %Sidelength=0.6;                       %���߳�
% 
% 
% position1=[-1,-0.5,0];
% f1=f;omiga1=f1*2*pi;wavek1=omiga1/speed;
% sourceInfo1=[position1,wavek1,omiga1,row,a,v0];
% 
% position2=[0,0.5,0];
% 
% f2=f+1000;omiga2=f2*2*pi;wavek2=omiga2/speed;
% sourceInfo2=[position2,wavek2,omiga2,row,a,v0];
% 
% NFFT=256;NBAND=NFFT/2.56+1;
% fs=16384;dt=1/fs;df=fs/NFFT;
% t=dt*[0:NFFT-1];fd=df*[0:NFFT-1];Nt=length(t);
% flow=f2-20; fhigh=f2+20;
% Nf1=round(flow/df)+1;Nf2=round(fhigh/df)+1;
% Nf=round(f/df)+1;
% 
% % ������grid type��������������������������������
% aperture size=0.8m,num of microphones is 35
% xh2=-0.3:0.1:0.3;
% yh2=-0.4:0.2:0.4;
% xh=xh2'*ones(1,length(yh2));xh=reshape(xh,1,[]);
% yh=yh2'*ones(1,length(xh2));yh=reshape(yh',1,[]);
% zh=zH*ones(1,length(xh));Nxh=length(xh);
% 
%% ������circle type������������������������������
% aperture size=0.2m,num of microphones is 360
theta=[0:359];
r=0.2;
xh=r*cos(theta*pi/180);yh=r*sin(theta*pi/180);Nxh=length(xh);
zh=zH*ones(1,length(xh));
% 
% %% ��������BK type����������������������������
% %aperture size=0.8m,num of microphones is 36
% % xh=[-0.001 -0.11 -0.218 -0.326 -0.262 -0.227 -0.111 -0.007 -0.023 -0.065 -0.005 0.09 0.044 -0.093 0.131 0.227 0.264 0.324 0.194 0.115 0.192 0.245 0.28 0.154 -0.043 0.017 0.087 0.098 -0.002 -0.907 -0.083 -0.108 -0.169 -0.179 -0.275 -0.244];
% % yh=[0.05,0.044,0.025,0.074,0.154,0.243,0.16,0.128,0.209,0.245,0.321,0.23,-0.025,0.074,0.176,0.245,0.15,0.075,0.016,0.058,-0.085,-0.066,-0.156,-0.193,-0.026,-0.117,-0.201,-0.319,-0.304,-0.319,-0.176,-0.07,-0.124,-0.179,-0.164,-0.037];
% % Nxh=length(xh);
% % zh=zH*ones(1,length(xh));
% 
% 
% 
%% ----------�ؽ���------------------------------------------
xr2=[-2:0.1:2];yr2=[-2:0.1:2];xrN=length(xr2);
xr=xr2'*ones(1,length(yr2));xr=reshape(xr,1,[]);
yr=yr2'*ones(1,length(xr2));yr=reshape(yr',1,[]);
zr=zR*ones(1,length(xr));Nxr=length(xr);
xrM=Nxr/xrN;    %--xrNΪ�еĳ��ȣ�xrMΪ�еĳ���

%% plot the microphone and sources topology
figure;
plot3(xh,yh,zh,'.b');
hold on;
plot3(xr,yr,zr,'Or','LineWidth',1.4);
plot3(realSrcCoord(1,1),realSrcCoord(1,2),realSrcCoord(1,3), '.r','markersize',30);
plot3(realSrcCoord(2,1),realSrcCoord(2,2),realSrcCoord(2,3), '.g','markersize',30);
plot3(realSrcCoord(3,1),realSrcCoord(3,2),realSrcCoord(3,3), '.m','markersize',30);



% 
% %% --��������------------------------------------------
% PH=BallSource2(sourceInfo1,t,xh,yh,zh,sourceInfo2);
% PH=PH.';
% % PH=BallSource2(sourceInfo1,t,xh,yh,zh);PH=PH.';
% 
% 
%% ---------�������þ�����ʽƵ�򷽷�����---------------

        R=[xr' yr' zr'];
        H=[xh' yh' zh'];
        xrr=reshape(xr,xrM,xrN);
        yrr=reshape(yr,xrM,xrN);
        Dist=sqrt(R.^2*ones(size(H'))+ones(size(R))*(H').^2-2*R*H');
        w=2*pi*freq';

      Nf1=  round(rotor_speed/60*29*1/df)-20;
      Nf2=  round(rotor_speed/60*29*1/df)+20;
        %ww����ɨ������,cormatrix������ؾ���       
        for n=1:Nxr
            rm=Dist(n,:);
            ww=exp(-1i*w.'*rm/speed)/Nxh;
%             for nf=1:NFFT         
%                cormatrix=H_fft(nf,:).'*conj(H_fft(nf,:));
% %                cutdiag=triu(cormatrix,1)+tril(cormatrix,-1);%ȥ�Խ�
%                temppow(nf)=conj(ww(nf,:))*cormatrix'*ww(nf,:).';%����㹫ʽ��ͬ�ĵط��ǣ�cormatrix����Ҫ����ת��    
%             end
             for nf=  Nf1:Nf2      %Nf1:Nf2 
               cormatrix=CC1(nf,:).'*conj(CC1(nf,:));
%                cutdiag=triu(cormatrix,1)+tril(cormatrix,1);%ȥ�Խ�
               temppow(nf)=conj(ww(nf,:))*cormatrix'*ww(nf,:).';%cormatrixδ����ת��    
            end           
            
            pow(n)=sum(temppow);%ȥ�ԽǵĹ���
            
        end
pow1=reshape(pow,xrM,xrN);
figure
contourf(xrr,yrr,abs(pow1));

