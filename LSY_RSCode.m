function encoded_msg = LSY_RSCode(dataIn_crc, nn, kk)
% RS编码函数
% 输入参数:
%   dataIn_crc: 输入的二进制数据
%   nn: 编码后码字长度(信息段+监督段)
%   kk: 信息段长度
% 输出参数:
%   encoded_msg: 编码后的二进制数据

% 初始化编码后的消息
encoded_msg = [];

% 计算每组处理的比特数
m = log2(nn + 1); % RS码中每个符号的比特数
bits_per_group = m * kk;

% 将输入数据分组处理
num_groups = floor(length(dataIn_crc) / bits_per_group);

for i = 1:num_groups
    % 提取当前组的数据
    start_idx = (i-1) * bits_per_group + 1;
    end_idx = i * bits_per_group;
    group_data = dataIn_crc(start_idx:end_idx);
    
    % 将二进制数据转换为GF(2^m)域中的符号
    gf_symbols = zeros(1, kk);
    for j = 1:kk
        bit_start = (j-1) * m + 1;
        bit_end = j * m;
        if bit_end <= length(group_data)
            symbol_bits = group_data(bit_start:bit_end);
            gf_symbols(j) = bi2de(symbol_bits, 'left-msb');
        end
    end
    
    % 创建RS编码器并编码
    rs_encoder = comm.RSEncoder('CodewordLength', nn, 'MessageLength', kk);
    encoded_symbols = rs_encoder(gf_symbols');
    
    % 将编码后的符号转换回二进制
    encoded_bits = [];
    for j = 1:nn
        symbol_dec = encoded_symbols(j);
        symbol_bits = de2bi(symbol_dec, m, 'left-msb');
        encoded_bits = [encoded_bits, symbol_bits];
    end
    
    % 添加到编码后的消息
    encoded_msg = [encoded_msg, encoded_bits];
end

% 处理剩余的比特
if length(dataIn_crc) > num_groups * bits_per_group
    remaining_bits = dataIn_crc(num_groups * bits_per_group + 1:end);
    encoded_msg = [encoded_msg, remaining_bits];
end

end

