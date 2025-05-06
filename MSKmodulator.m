function [biNRZ, txSig] = MSKmodulator(samp, TX_BIT_MAT)
% MSK调制器函数
% 输入参数:
%   samp: 过采样倍数
%   TX_BIT_MAT: 输入的二进制数据矩阵，每行对应一个跳频点的数据
% 输出参数:
%   biNRZ: 双极性NRZ信号
%   txSig: MSK调制后的信号

% 获取矩阵维度
[num_hops, bits_per_hop] = size(TX_BIT_MAT);

% 初始化输出
biNRZ = zeros(num_hops, bits_per_hop);
txSig = zeros(num_hops, samp * bits_per_hop);

% MSK参数
h = 0.5;  % 调制指数，MSK固定为0.5

for hop = 1:num_hops
    % 获取当前跳频点的比特
    bits = TX_BIT_MAT(hop, :);
    
    % 转换为双极性NRZ信号 (0->-1, 1->1)
    biNRZ(hop, :) = 2 * bits - 1;
    
    % 对每个比特进行过采样
    biNRZ_oversamp = zeros(1, samp * bits_per_hop);
    for i = 1:bits_per_hop
        biNRZ_oversamp((i-1)*samp+1:i*samp) = biNRZ(hop, i);
    end
    
    % MSK调制
    % 初始化I和Q分量
    I_phase = zeros(1, samp * bits_per_hop);
    Q_phase = zeros(1, samp * bits_per_hop);
    
    % 生成I和Q分量的相位
    for i = 1:bits_per_hop
        bit_idx = (i-1)*samp+1:i*samp;
        t = linspace(0, 1, samp);
        
        if mod(i, 2) == 1  % 奇数位影响I分量
            I_phase(bit_idx) = biNRZ(hop, i) * cos(pi*t/2);
        else  % 偶数位影响Q分量
            Q_phase(bit_idx) = biNRZ(hop, i) * sin(pi*t/2);
        end
    end
    
    % 合成MSK信号
    carrier_I = cos(2*pi*0.25*(0:samp*bits_per_hop-1)/samp);
    carrier_Q = sin(2*pi*0.25*(0:samp*bits_per_hop-1)/samp);
    
    I_signal = I_phase .* carrier_I;
    Q_signal = Q_phase .* carrier_Q;
    
    % 复信号表示
    txSig(hop, :) = I_signal + 1j * Q_signal;
end

end

