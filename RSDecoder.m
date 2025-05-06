function dataOut_dc = RSDecoder(deinterwine_data, nn, kk)
% RS解码函数
% 输入参数:
%   deinterwine_data: 输入的二进制数据
%   nn: 编码后码字长度(信息段+监督段)
%   kk: 信息段长度
% 输出参数:
%   dataOut_dc: 解码后的二进制数据

% 初始化解码后的消息
dataOut_dc = [];

% 计算每组处理的比特数
m = log2(nn + 1); % RS码中每个符号的比特数
bits_per_codeword = m * nn;

% 将输入数据分组处理
num_codewords = floor(length(deinterwine_data) / bits_per_codeword);

for i = 1:num_codewords
    % 提取当前码字的数据
    start_idx = (i-1) * bits_per_codeword + 1;
    end_idx = i * bits_per_codeword;
    codeword_bits = deinterwine_data(start_idx:end_idx);
    
    % 将二进制数据转换为GF(2^m)域中的符号
    gf_symbols = zeros(1, nn);
    for j = 1:nn
        bit_start = (j-1) * m + 1;
        bit_end = j * m;
        if bit_end <= length(codeword_bits)
            symbol_bits = codeword_bits(bit_start:bit_end);
            gf_symbols(j) = bi2de(symbol_bits, 'left-msb');
        end
    end
    
    % 创建RS解码器并解码
    rs_decoder = comm.RSDecoder('CodewordLength', nn, 'MessageLength', kk);
    try
        % 直接使用数值，而不是gf对象
        decoded_symbols = rs_decoder(gf_symbols');
        
        % 将解码后的符号转换回二进制
        decoded_bits = [];
        for j = 1:kk
            symbol_dec = decoded_symbols(j);
            symbol_bits = de2bi(symbol_dec, m, 'left-msb');
            decoded_bits = [decoded_bits, symbol_bits];
        end
    catch
        % 如果解码失败，保留原始信息部分
        decoded_bits = codeword_bits(1:m*kk);
    end
    
    % 添加到解码后的消息
    dataOut_dc = [dataOut_dc, decoded_bits];
end

% 处理剩余的比特
if length(deinterwine_data) > num_codewords * bits_per_codeword
    remaining_bits = deinterwine_data(num_codewords * bits_per_codeword + 1:end);
    dataOut_dc = [dataOut_dc, remaining_bits];
end

end
