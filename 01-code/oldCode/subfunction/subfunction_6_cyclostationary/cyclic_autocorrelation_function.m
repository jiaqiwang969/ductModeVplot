function R = cyclic_autocorrelation_function(x,Nfft,max_tau)
%���룺x��һάʱ������
%      Nfft�����ڿ��ٸ���Ҷ�任�����ݵ���������ѭ��Ƶ�ʵ����ݵ���
%      max_tau�����ʱ�ӵ���
%�����ʱ������x��ѭ������غ�������R
%      Rά��: [2* max_tau+1, Nfft]

[rows,cols] = size(x);
if rows > cols
  x = x';
end
x = x - mean(x);

n = floor((length(x)-2*max_tau-1)/Nfft);  % ȷ����ƽ���Ĵ���
r = zeros(2*max_tau+1,Nfft);           % ʱ������غ���
temp = zeros(Nfft,n);
t = (1:Nfft*n); 
for k = -max_tau:max_tau   
   temp(:) = conj(x(t+max_tau)).*x(t+k+max_tau);  %˳������������
   if n==1
        r(k+1+max_tau,:) = temp';
    else
        r(k+1+max_tau,:) = mean(temp');   %��ÿ��������ƽ��
   end 
end

R = zeros(2*max_tau+1,Nfft);
for k = -max_tau:max_tau
R(k+1+max_tau,:) = exp(-j*pi*((0:Nfft-1)/Nfft)*k).*fft(r(k+1+max_tau,:))/Nfft;  
%�Գ���ʽ��ѭ������غ���
end

R = fftshift(R,2);  %��ѭ��Ƶ����0Ϊ���ĵ�ԳƷ���
