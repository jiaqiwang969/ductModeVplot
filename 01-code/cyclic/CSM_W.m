function Mat = CSM_W(y,x,nfft,Noverlap,Window)     
% Mat = CSM_W(y,x,nfft,Noverlap,Window)
% �ź� x �� �ź� y ֮��Ļ��׾������(����)
%
% ����:
% y, x = �źž���, ����ÿһ����һ��ͨ�����ź�
% Window = ���ڳ���
% nfft = ÿ��������FFT��ĸ��� �Ƽ���nfft = 2*NWind, �����ڳ���NWind������
% Noverlap = �����ص� �Ƽ���Noverlap = 1/2*Nwind (������)�� Noverlap = 2/3*Nwind (������)
% CSM_W ���� 'CPS_W' ����.
%
% ���: 
% S ��һ���ṹ���͵�����:
%                   Mat.S = ���׾���, ����ά��Ϊ (x��)x(y��)x(Nwind/2+1)
%                   Mat.f  = Ƶ������
%                   Mat.K = ����block�ĸ���
%                   Mat.Var_Reduc = �����С����


Kx = size(x,2);
Ky = size(y,2);

Nw = length(Window);
if Nw == 1
    Nw = Window;
end
Mat.S = zeros(Ky,Kx,Nw/2+1);
count = 1;

same = 0;
if Kx == Ky
    if all(all(y==x))
        same = 1;
    end
end
if same == 1
    for i = 1:Kx
        Sij = CPS_W(y(:,i),x(:,i),nfft,Noverlap,Window);
        Mat.S(i,i,:) = real(Sij.S(1:Nw/2+1));
        for j = 1:i-1
            waitbar(count/Kx/(Kx-1)*2),count = count + 1;
            Sij = CPS_W(y(:,i),x(:,j),nfft,Noverlap,Window);
            Mat.S(i,j,:) = Sij.S(1:Nw/2+1);
            Mat.S(j,i,:) = conj(Sij.S(1:Nw/2+1));
        end
    end
else
    for i = 1:Ky
        for j = 1:Kx
            waitbar(count/Kx/Ky),count = count + 1;
            Sij = CPS_W(y(:,i),x(:,j),nfft,Noverlap,Window);
            Mat.S(i,j,:) = Sij.S(1:Nw/2+1);
        end
    end
end

Mat.Var_Reduc = Sij.Var_Reduc; 
Mat.K = Sij.K;
Mat.f = Sij.f; 


function Spec = CPS_W(y,x,nfft,Noverlap,Window)
% Spec = CPS_W(y,x,Noverlap,Window)
% �ź� x �� y ֮��Ļ��׹��� (����)
% ����:
% y, x = �����ź�
% Window = ���ڳ���
% nfft = ÿ��������FFT��ĸ��� �Ƽ���nfft = 2*NWind, �����ڳ���NWind������
% Noverlap = �����ص� �Ƽ���Noverlap = 1/2*Nwind (������)�� Noverlap = 2/3*Nwind (������)
%
% ���:
% Spec �ǽṹ���͵��������:
%                   Spec.S = ��������
%                   Spec.f =  Ƶ��(����)
%                   Spec.K = ����block�ĸ���
%                   Spec.Var_Reduc = �����С����

if length(Window) == 1
    Window = hanning(Window);
end
Window = Window(:);
n = length(x);          % ���ݵ����
nwind = length(Window); % ���ڳ���

% ������
if nwind <= Noverlap,error('Window length must be > Noverlap');end
if nfft < nwind,error('Window length must be <= nfft');end

y = y(:);
x = x(:);
K = fix((n-Noverlap)/(nwind-Noverlap));	% ���ڸ���

% ���� CPS ��
index = 1:nwind;
f = (0:nfft-1)/nfft;
CPS = 0;

for i=1:K
    Yw = fft(Window.*y(index),nfft);		% Yw(f)
    Xw = fft(Window.*x(index),nfft);		% Xw(f-a)
    CPS = Yw.*conj(Xw) + CPS;
    index = index + (nwind - Noverlap);
end

% ��һ��
KMU = K*norm(Window)^2;	% ��һ������, ��֤������ƫ
CPS = CPS/KMU;

% ���㷽���С����
Window = Window(:)/norm(Window);
Delta = nwind - Noverlap;
R2w = xcorr(Window);
k = nwind+Delta:Delta:min(2*nwind-1,nwind+Delta*(K-1));
if length(k) >1
    Var_Reduc = R2w(nwind)^2/K + 2/K*(1-(1:length(k))/K)*(R2w(k).^2);
else
    Var_Reduc = R2w(nwind)^2/K;
end

% �������
Spec.S = CPS;
Spec.f = f;
Spec.K = K;
Spec.Var_Reduc = Var_Reduc;





