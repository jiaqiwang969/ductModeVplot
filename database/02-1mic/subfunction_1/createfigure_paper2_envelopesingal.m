function createfigure_paper2_envelopesingal(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  x 数据的矢量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 12-Jul-2019 10:32:24 自动生成

% 创建 figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% 创建 axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.57152855993564 0.177744227353464]);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1);

% 取消以下行的注释以保留坐标轴的 X 范围
 xlim(axes1,[10.7247422903382 200]);
% 取消以下行的注释以保留坐标轴的 Y 范围
ylim(axes1,[-8.26559906251832e-06 0.15]);
% 取消以下行的注释以保留坐标轴的 Z 范围
% zlim(axes1,[-1 1]);
box(axes1,'on');
grid(axes1,'on');
% 设置其余坐标轴属性
set(axes1,'GridLineStyle','--');
