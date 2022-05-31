function  [GAMMA_residual,freq_residual]=CC(y_residual)
%% 方法02:合成孔径算法
zH = 0;             % 测试距离
zRef=zH-0.08;       % 参考传声器到测试面的轴向距离
NumMic = 12;        % 传声器的数量
NumSM  = 30;        % 测量的次数
Radius = 0.2;       % 管道半径
Area = pi*Radius^2; % 管口面积
Fs = 102400;        % 采样频率
time = 5;           % 采样时间

L_signal = Fs*time;             %信号长度
L_seg = round(L_signal/100);    %确定对信号处理划窗长度
Wind =  hamming(L_seg);         %确定对数据进行汉宁窗处理
Noverlap = round(L_seg/2);      %确定信号划窗重叠率
Nfft = 2^(ceil(log2(L_seg))+1); %确定分析频率
block=16;%分段,生成多个block，每个block 为round

%% 方法01:CC算法
Ind = [1:NumSM];   %设定循环次数
Num_file = Ind ;
for i_file =Num_file
    Tdata_residual{i_file}=y_residual(:,13*i_file-12:13*i_file);
    [temp_ref_residual,freq_residual] = cpsd(Tdata_residual{i_file}(:,13),Tdata_residual{i_file}(:,13),Wind,Noverlap,Nfft,Fs);
    temp_ref_residual = sqrt(temp_ref_residual);

    for k=1:NumMic
        [temp_residual,freq_residual] = cpsd(Tdata_residual{i_file}(:,k),Tdata_residual{i_file}(:,13),Wind,Noverlap,Nfft,Fs);
        CC_residual(:, i_file+NumSM*(k-1)) = temp_residual./temp_ref_residual;  %360个测点的数据矩阵
    end
end

NumMic_enlarge=NumSM*NumMic;                                                                             %传声器可测量的模态总数
mode=[-60:60];
for k=1:length(mode)
    a_mf_residual(:,k)=1/NumMic_enlarge*CC_residual(:,1:NumMic_enlarge)*exp(mode(k)*1i*2*pi*(0:NumMic_enlarge-1)/NumMic_enlarge).';    %求解模态系数

end
GAMMA_residual =abs(a_mf_residual);

end
