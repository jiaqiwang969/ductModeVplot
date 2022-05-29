function createfigure_paper2_envelopesingal(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  x ���ݵ�ʸ��
%  YMATRIX1:  y ���ݵľ���

%  �� MATLAB �� 12-Jul-2019 10:32:24 �Զ�����

% ���� figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% ���� axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.57152855993564 0.177744227353464]);
hold(axes1,'on');

% ʹ�� plot �ľ������봴������
plot1 = plot(X1,YMatrix1);

% ȡ�������е�ע���Ա���������� X ��Χ
 xlim(axes1,[10.7247422903382 200]);
% ȡ�������е�ע���Ա���������� Y ��Χ
ylim(axes1,[-8.26559906251832e-06 0.15]);
% ȡ�������е�ע���Ա���������� Z ��Χ
% zlim(axes1,[-1 1]);
box(axes1,'on');
grid(axes1,'on');
% ������������������
set(axes1,'GridLineStyle','--');
