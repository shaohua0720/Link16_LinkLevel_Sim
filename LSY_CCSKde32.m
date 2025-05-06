function ccsk_dec = LSY_CCSKde32(dmdBit, ccskcode)
% CCSK32解扩函数
% 输入参数:
%   dmdBit: 输入的扩频后二进制数据
%   ccskcode: CCSK序列(32位)
% 输出参数:
%   ccsk_dec: 解扩后的二进制数据

% CCSK序列长度
n_ccsk = length(ccskcode);

% 计算输入数据中包含的符号数
num_symbols = floor(length(dmdBit) / n_ccsk);

% 初始化解扩后的数据
ccsk_dec = zeros(1, num_symbols);

% 生成循环移位后的CCSK序列（用于比较）
shifted_code = [ccskcode(end); ccskcode(1:end-1)];

% 对每个CCSK符号进行解扩
for i = 1:num_symbols
    % 提取当前CCSK符号
    start_idx = (i-1) * n_ccsk + 1;
    end_idx = i * n_ccsk;
    current_symbol = dmdBit(start_idx:end_idx);
    
    % 计算与原始CCSK序列的相关性
    corr_original = sum(current_symbol == ccskcode');
    
    % 计算与移位CCSK序列的相关性
    corr_shifted = sum(current_symbol == shifted_code');
    
    % 根据相关性判断原始比特
    if corr_shifted > corr_original
        ccsk_dec(i) = 1;
    else
        ccsk_dec(i) = 0;
    end
end

% 处理剩余的比特
remaining_bits = length(dmdBit) - num_symbols * n_ccsk;
if remaining_bits > 0
    ccsk_dec = [ccsk_dec, dmdBit(num_symbols * n_ccsk + 1:end)];
end

end