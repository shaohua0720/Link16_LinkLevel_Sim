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

% 将信号重塑为每比特samp个样本
rxSamples = reshape(rxSig(1:bitCount*samp), samp, bitCount);

% 分离I和Q分量
I_samples = real(rxSamples);
Q_samples = imag(rxSamples);

% 生成本地载波
t = (0:samp-1) / samp;
carrier_I = cos(2*pi*0.25*t);
carrier_Q = sin(2*pi*0.25*t);

% 解调
for i = 1:bitCount
    % 获取当前比特的样本
    I_bit = I_samples(:, i)';
    Q_bit = Q_samples(:, i)';
    
    % 解调I和Q分量
    I_demod = I_bit .* carrier_I;
    Q_demod = Q_bit .* carrier_Q;
    
    % 计算I和Q分量的相关性
    I_corr = sum(I_demod .* cos(pi*t/2));
    Q_corr = sum(Q_demod .* sin(pi*t/2));
    
    % 根据比特位置(奇/偶)决定使用I或Q分量
    if mod(i, 2) == 1  % 奇数位
        demodBits(i) = I_corr > 0;
    else  % 偶数位
        demodBits(i) = Q_corr > 0;
    end
end

end