function rcvBBmat = FHdemodulator(rcvNoisyMat, txFHtable, fs)
% 跳频解调器函数
% 输入参数:
%   rcvNoisyMat: 接收到的带噪声的跳频信号矩阵
%   txFHtable: 跳频频点表
%   fs: 采样率
% 输出参数:
%   rcvBBmat: 解跳频后的基带信号矩阵

% 获取跳频点数和每跳的采样点数
[HOP_NUM, samplesPerHop] = size(rcvNoisyMat);

% 初始化输出矩阵
rcvBBmat = zeros(HOP_NUM, samplesPerHop);

% 对每个跳频点进行处理
for hop = 1:HOP_NUM
    % 获取当前跳频点的频率
    freq = txFHtable(hop);
    
    % 获取当前跳频点的接收信号
    received_signal = rcvNoisyMat(hop, :);
    
    % 生成时间向量
    t = (0:samplesPerHop-1) / fs;
    
    % 下变频到基带
    carrier = exp(-1j * 2 * pi * freq * t);
    
    % 跳频解调
    rcvBBmat(hop, :) = received_signal .* carrier;
    
    % 可选：添加低通滤波器以去除高频分量
    % rcvBBmat(hop, :) = lowpass(rcvBBmat(hop, :), cutoff_freq, fs);
end

end