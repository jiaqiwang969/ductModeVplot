function [up,down] = envelope_universal(x,y,interpMethod)
%https://zhuanlan.zhihu.com/p/21335419
%   ------------------------------------------
%   [up,down] = envelope(x,y,interpMethod)
%
%   ����:
%    x               ���ݺ�����
%    y               ����������
%    interpMethod    ��ֵ����
%
%   ���:
%    up      �ϰ��磬��xͬ���ĳ���
%    down    �°��磬��xͬ���ĳ���

if length(x) ~= length(y)
   error('Two input data should have the same length��');
end

if (nargin < 2)||(nargin > 3),
   error('Please see help for INPUT DATA.');
elseif (nargin == 2)
   interpMethod = 'linear';
end

%-----------------------------------------------  
% �����ֵ�����Ӧ����
extrMaxValue = y(find(dif(sign(dif(y)))==-2)+1); 
extrMaxIndex =   find(dif(sign(dif(y)))==-2)+1;  

%-----------------------------------------------
% ����С��ֵ�����Ӧ����
extrMinValue = y(find(dif(sign(dif(y)))== +2)+1); 
extrMinIndex =   find(dif(sign(dif(y)))== +2)+1;  

up = extrMaxValue;
up_x = x(extrMaxIndex);

down = extrMinValue;
down_x = x(extrMinIndex);

% -----------------------------------------
% ���°����ֵ
up = interp1(up_x,up,x,interpMethod); 
down = interp1(down_x,down,x,interpMethod); 

% --------------------------
% ���Ǳ߽����������
up(1) = y(1);       % �ױ߽�
up(end) = y(end);   % β�߽�
down(1) = y(1);     % �ױ߽�
down(end) = y(end); % β�߽�
N = length(x);
for i = 1:N
    if isinf(up(i))==1
        up(i) = y(i);
    end
    if isinf(down(i))==1
        down(i) = y(i);
    end
end


function dy = dif(y)
dy = zeros(size(y));
dy(1:end-1) = y(2:end) - y(1:end-1);
dy(end) = (y(end-2) + 2*y(end-1)- 3*y(end))/6;