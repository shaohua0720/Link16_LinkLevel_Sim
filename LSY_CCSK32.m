function ccsk_msg = LSY_CCSK32(interwine_msg, ccskcode)
% CCSK32扩频函数
% 输入参数:
%   interwine_msg: 输入的二进制数据
%   ccskcode: CCSK序列(32位)
% 输出参数:
%   ccsk_msg: 扩频后的二进制数据

% 初始化扩频后的消息
ccsk_msg = [];

% CCSK序列长度
n_ccsk = length(ccskcode);

% 处理每个输入比特
for i = 1:length(interwine_msg)
    if interwine_msg(i) == 0
        % 如果输入比特为0，直接使用CCSK序列
        spread_bits = ccskcode';
    else
        % 如果输入比特为1，使用CCSK序列的循环移位版本
        % Link16使用的是循环右移1位
        shifted_code = [ccskcode(end); ccskcode(1:end-1)];
        spread_bits = shifted_code';
    end
    
    % 添加到扩频后的消息
    ccsk_msg = [ccsk_msg, spread_bits];
end

end