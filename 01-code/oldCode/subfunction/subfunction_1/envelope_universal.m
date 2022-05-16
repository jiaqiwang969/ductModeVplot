function [up,down] = envelope_universal(x,y,interpMethod)
%https://zhuanlan.zhihu.com/p/21335419
%   ------------------------------------------
%   [up,down] = envelope(x,y,interpMethod)
%
%   输入:
%    x               数据横坐标
%    y               数据纵坐标
%    interpMethod    插值方法
%
%   输出:
%    up      上包络，跟x同样的长度
%    down    下包络，跟x同样的长度

if length(x) ~= length(y)
   error('Two input data should have the same length！');
end

if (nargin < 2)||(nargin > 3),
   error('Please see help for INPUT DATA.');
elseif (nargin == 2)
   interpMethod = 'linear';
end

%-----------------------------------------------  
% 众最大极值及其对应坐标
extrMaxValue = y(find(dif(sign(dif(y)))==-2)+1); 
extrMaxIndex =   find(dif(sign(dif(y)))==-2)+1;  

%-----------------------------------------------
% 众最小极值及其对应坐标
extrMinValue = y(find(dif(sign(dif(y)))== +2)+1); 
extrMinIndex =   find(dif(sign(dif(y)))== +2)+1;  

up = extrMaxValue;
up_x = x(extrMaxIndex);

down = extrMinValue;
down_x = x(extrMinIndex);

% -----------------------------------------
% 上下包络插值
up = interp1(up_x,up,x,interpMethod); 
down = interp1(down_x,down,x,interpMethod); 

% --------------------------
% 考虑边界与特殊情况
up(1) = y(1);       % 首边界
up(end) = y(end);   % 尾边界
down(1) = y(1);     % 首边界
down(end) = y(end); % 尾边界
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