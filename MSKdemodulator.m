function [demodBits] = MSKdemodulator(samp, rxSig)
% MSK解调器函数
% 输入参数:
%   samp: 过采样倍数
%   rxSig: 接收到的MSK调制信号
% 输出参数:
%   demodBits: 解调后的二进制数据

% 确保rxSig是行向量
rxSig = rxSig(:)';

% 计算信号长度和比特数
sigLen = length(rxSig);
bitCount = floor(sigLen / samp);

% 初始化输出
demodBits = zeros(1, bitCount);

% 提取I和Q分量
I_signal = real(rxSig);
Q_signal = imag(rxSig);

% 生成本地载波
t_full = (0:sigLen-1) / samp;
carrier_I = cos(2*pi*0.25*t_full);
carrier_Q = sin(2*pi*0.25*t_full);

% 下变频
I_baseband = I_signal .* carrier_I;
Q_baseband = Q_signal .* carrier_Q;

% 对每个比特进行处理
for i = 1:bitCount
    % 计算当前比特的样本索引
    bit_idx = (i-1)*samp+1:i*samp;
    t_bit = (0:samp-1) / samp;
    
    if mod(i, 2) == 1  % 奇数位 - 影响I分量
        % 使用余弦波形进行相关
        I_bit = I_baseband(bit_idx);
        corr_value = sum(I_bit .* cos(pi*t_bit/2));
        demodBits(i) = corr_value > 0;
    else  % 偶数位 - 影响Q分量
        % 使用正弦波形进行相关
        Q_bit = Q_baseband(bit_idx);
        corr_value = sum(Q_bit .* sin(pi*t_bit/2));
        % 修改：反转偶数位的判断条件
        demodBits(i) = corr_value < 0;  % 注意这里改为 < 0
    end
end

end


