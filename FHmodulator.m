function txFHmodulatedMat = FHmodulator(msgModMatrix, txFHtable, fs)
% 跳频调制器函数
% 输入参数:
%   msgModMatrix: 调制后的信号矩阵，每行对应一个跳频点的数据
%   txFHtable: 跳频频点表
%   fs: 采样率
% 输出参数:
%   txFHmodulatedMat: 跳频调制后的信号矩阵

% 获取跳频点数和每跳的采样点数
[HOP_NUM, samplesPerHop] = size(msgModMatrix);

% 初始化输出矩阵
txFHmodulatedMat = zeros(HOP_NUM, samplesPerHop);

% 对每个跳频点进行处理
for hop = 1:HOP_NUM
    % 获取当前跳频点的频率
    freq = txFHtable(hop);
    
    % 获取当前跳频点的基带信号
    baseband_signal = msgModMatrix(hop, :);
    
    % 生成时间向量
    t = (0:samplesPerHop-1) / fs;
    
    % 上变频到指定频点
    carrier = exp(1j * 2 * pi * freq * t);
    
    % 跳频调制
    txFHmodulatedMat(hop, :) = baseband_signal .* carrier;
end

end