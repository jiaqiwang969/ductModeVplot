function [y,g] = Filt_STFT(x,G)
% [y,g] = Filt_STFT(x,G)
%
% Filters signals x with the frequency response matrix G to give signals y.
% G has dimensions (nb of signals y)*(nb of signals x)*(nb of freq. bins).
% g is the returned time-filter corresponding to G.
%


[L,mx] = size(x);
if mx > L,x = x.';[L,mx]=size(x);end
N = size(G);
L_N = length(N);
IL = find(N==max(N));
if L_N == 3
   Imx = find(N==mx);
   Imy = setdiff(1:L_N,[Imx,IL]);
   G = permute(G,[Imy,Imx,IL]);
elseif mx == 1
   Imy = setdiff(1:L_N,IL);
   G = permute(G,[Imy,IL]);
else
   Imx = setdiff(1:L_N,IL);
   G = permute(G,[Imx,IL]);
end
N = size(G);
g = zeros(N);

if L_N == 3
   if N(2) ~= mx,error('Dimension mismatch between x and G !!'),end
end  

N2 = floor(N(end)/2);
if L_N == 3
   y = zeros(L+N2,N(1));
	for i=1:N(1)
        for j=1:N(2)
            g(i,j,:) = fftshift(real(ifft(G(i,j,:))));
            y(:,i) = y(:,i) + fftfilt(squeeze(g(i,j,:)),[x(:,j);zeros(N2,1)]);
        end
	end
elseif mx == 1
   y = zeros(L+N2,N(1));
   for i=1:N(1)
      g(i,:) = fftshift(real(ifft(G(i,:))));
      y(:,i) = fftfilt(g(i,:),[x;zeros(N2,1)]);
   end
else
   y = 0;
   for j=1:N(1)
      g(j,:) = fftshift(real(ifft(G(j,:))));
      y = y + fftfilt(g(j,:),[x(:,j);zeros(N2,1)]);
   end
end
y = y(N2+1:end,:);
