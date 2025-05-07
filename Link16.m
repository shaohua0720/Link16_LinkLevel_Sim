%% Link16通信系统仿真
% 版本: v1.0.0
% 跳频通信系统误比特率统计分析
% 处理流程: 调制-->跳频-->信道-->解跳-->解调-->误码分析
tic
clc; clear; close all;

%% 1. 系统参数设置
% ===== 1.1 系统功能开关 =====
ENABLE_FH = 0;                  % 0: 不跳频传输 1: 跳频传输
MODULATION_TYPE = 2;            % 0: BPSK, 1: QPSK, 2: 16QAM, 3: 64QAM, 4: MSK
ENABLE_CRC = 0;                 % 0: 不使用CRC编码, 1: 使用CRC编码
CODING_TYPE = 1;                % 0: 不编码, 1: RS(31,15)编码, 2: (2,1,7)卷积编码
ENABLE_INTERLEAVING = 1;        % 0: 不使用交织, 1: 使用交织
ENABLE_CCSK = 1;                % 0: 不使用CCSK扩频, 1: 使用CCSK扩频

% ===== 1.2 基本系统参数 =====
bitRate = 5e4;                  % 比特率: 50Kb/s
samplesPerBit = 20;             % 过采样倍数
samplingRate = samplesPerBit * bitRate;  % 采样率: 1MHz

% ===== 1.3 跳频参数 =====
hoppingRate = 1000;             % 跳频速率: 1000跳/秒
bitsPerHop = bitRate / hoppingRate;  % 每跳比特数: 50比特/跳
frequencyBandwidth = 5e6;       % 跳频带宽: 5MHz
frequencyPointCount = floor(frequencyBandwidth / (bitRate*4));  % 跳频频点数: 25个频点
frequencyInterval = frequencyBandwidth / frequencyPointCount;  % 频点间隔: 200kHz
frequencySequence = ((0:frequencyPointCount-1) - floor(frequencyPointCount/2)) * frequencyInterval;  % 跳频频点序列: 相对于中心频率的偏移量
carrierFrequency = 300e6;       % 载波中心频率: 300MHz
carrierSequence = carrierFrequency + frequencySequence;  % 跳频载波频率序列: 绝对频率值

% ===== 1.4 调制参数 =====
switch MODULATION_TYPE
    case 0  % BPSK
        bitsPerSymbol = 1;      % BPSK: 1比特/调制符号
        figureLabel_Modulation = 'BPSK';
    case 1  % QPSK
        bitsPerSymbol = 2;      % QPSK: 2比特/调制符号
        figureLabel_Modulation = 'QPSK';
    case 2  % 16QAM
        bitsPerSymbol = 4;      % 16QAM: 4比特/调制符号
        figureLabel_Modulation = '16QAM';
    case 3  % 64QAM
        bitsPerSymbol = 6;      % 64QAM: 6比特/调制符号
        figureLabel_Modulation = '64QAM';
    case 4  % MSK
        bitsPerSymbol = 1;      % MSK: 1比特/调制符号
        figureLabel_Modulation = 'MSK';
end
modulationOrder = 2^bitsPerSymbol;  % 调制阶数: 2^比特数

% ===== 1.5 CRC参数 =====
if ENABLE_CRC == 1
    crcBlockSize = 225;         % CRC分组大小: 每225比特信息添加校验
    figureLabel_CRC = 'CRC(237,225)';  % 225位信息 + 12位CRC校验
    polynomial = 'z12+1';       % CRC多项式: z^12+1
else
    crcBlockSize = 1;           % 不使用CRC时设为1以避免影响计算
    figureLabel_CRC = 'no CRC';
end

% ===== 1.6 编码参数 =====
rs_m = 5;                       % RS码符号位宽: 每个符号用5比特表示
rs_n = 2^rs_m - 1;              % RS码字长度: 31个符号(155比特)
rs_t = 8;                       % RS码纠错能力: 最多纠正8个符号错误
rs_k = rs_n - 2*rs_t;           % RS信息段长度: 15个符号(75比特)

switch CODING_TYPE
    case 0  % 不编码
        figureLabel_Coding = '无纠错编码';
    case 1  % RS(31,15)编码
        figureLabel_Coding = 'RS(31,15)编码';
    case 2  % (2,1,7)卷积编码
        figureLabel_Coding = '(2,1,7)卷积编码';
        constraintLength = 7;   % 卷积码约束长度: 7
        tracebackDepth = 42;    % Viterbi译码器回溯深度: 42 (约束长度的6倍)
        trellis = poly2trellis(constraintLength, [171, 133]);  % 卷积码生成多项式: 171和133 (八进制)
end

% ===== 1.7 交织参数 =====
if ENABLE_INTERLEAVING
    interleavingRows = 10;      % 交织矩阵行数: 10行
    interleavingCols = 100;     % 交织矩阵列数: 100列
    figureLabel_Interleaving = '使用交织';
else
    figureLabel_Interleaving = '不使用交织';
end

% ===== 1.8 CCSK参数 =====
if ENABLE_CCSK == 1
    % 优化的32位CCSK序列
    ccskCode = [1 0 1 1 1 0 1 0 0 0 1 1 1 1 0 1 0 0 1 0 0 0 0 0 0 1 1 0 0 1 1 0]';
    figureLabel_CCSK = 'CCSK';
else
    figureLabel_CCSK = '不使用CCSK';
end

% ===== 1.9 信息参数 =====
% 计算基本消息长度(考虑RS编码和调制要求)
messageLength = rs_m * rs_k;    % 基本消息长度: 75比特(RS编码的一个信息段)
frameCount = 1;                 % 传输帧数: 1帧
packetCount = 3;                % 每帧包数: 3包
bitsPerPacket = messageLength * bitsPerHop;  % 每包比特数: 75 × 50 = 3750比特
messageBitCount = bitsPerPacket * packetCount;  % 总消息比特数: 3750 × 3 = 11250比特
totalBitCount = messageBitCount;  % 总发送比特数: 11250比特

% ===== 1.10 跳频标签设置 =====
if ENABLE_FH == 1
    figureLabel_FH = '跳频传输';
else
    figureLabel_FH = '非跳频传输';
end

% ===== 1.11 创建MSK调制器和解调器对象 =====
if MODULATION_TYPE == 4  % MSK
    mskModulator = comm.MSKModulator('BitInput', true, ...
                                     'SamplesPerSymbol', samplesPerBit, ...
                                     'InitialPhaseOffset', 0);  % 设置初始相位为0

    mskDemodulator = comm.MSKDemodulator('BitOutput', true, ...
                                         'SamplesPerSymbol', samplesPerBit, ...
                                         'InitialPhaseOffset', 0);  % 使用相同的初始相位
end

% ===== 1.12 误比特分析参数 =====
snrValues = (-20:-10);            % Eb/No范围: -10dB到0dB，步长1dB
codeRate = 1/2;                 % 系统总码率: 1/2 (考虑编码和扩频)
berValues = zeros(1, length(snrValues));  % 误比特率结果数组
packetErrorRate = zeros(1, length(snrValues));  % 误包率结果数组
crcErrorCount = zeros(1, length(snrValues));  % CRC检测错误计数

%% 2. 发送端信号处理
% ===== 2.1 构造发送序列 =====
messageBits = randi([0, 1], 1, totalBitCount);  % 随机生成二进制消息
inputBits = messageBits;        % 构造发送比特序列

% ===== 2.2 CRC检错编码 =====
if ENABLE_CRC
    checksumPerFrame = length(inputBits)/crcBlockSize;  % CRC校验块数
    crcGenerator = comm.CRCGenerator(polynomial, 'ChecksumsPerFrame', checksumPerFrame);  % CRC编码器
    crcDetector = comm.CRCDetector(polynomial, 'ChecksumsPerFrame', checksumPerFrame);  % CRC检测器
    crcEncodedBits = crcGenerator(inputBits');  % CRC编码(输入为列向量)
    crcEncodedBits = crcEncodedBits';  % 转换回行向量
else
    crcEncodedBits = inputBits;  % 不使用CRC时直接传递
end

% ===== 2.3 信源纠错编码 =====
switch CODING_TYPE
    case 0  % 不编码
        codedBits = crcEncodedBits;
    case 1  % RS编码
        codewordLength = rs_n;  % RS码字长度: 31符号
        messageLength = rs_k;   % RS信息长度: 15符号
        codedBits = RSEncoder(crcEncodedBits, codewordLength, messageLength);  % RS编码
    case 2  % 卷积编码
        codedBits = convenc(crcEncodedBits, trellis);  % 卷积编码(码率1/2)
end

% ===== 2.4 交织 =====
if ENABLE_INTERLEAVING
    interleavingCount = length(codedBits)/(interleavingRows*interleavingCols);  % 交织块数
    interleavedBits = zeros(1, length(codedBits));

    for i = 1:interleavingCount
        startIdx = ((i-1)*(interleavingRows*interleavingCols))+1;
        endIdx = i*(interleavingRows*interleavingCols);
        tempData = codedBits(1, startIdx:endIdx);
        % 按列填充矩阵，按行读出实现交织
        interleavedData = matintrlv(tempData, interleavingRows, interleavingCols);
        interleavedBits(1:i*interleavingRows*interleavingCols) = ...
            horzcat(interleavedBits(1:(i-1)*interleavingRows*interleavingCols), interleavedData);
    end
else
    interleavedBits = codedBits;  % 不使用交织时直接传递
end

% ===== 2.5 CCSK扩频 =====
if ENABLE_CCSK == 1
    spreadBits = CCSKModulator(interleavedBits, ccskCode);  % CCSK扩频: 1比特→32比特
else
    spreadBits = interleavedBits;  % 不使用CCSK时直接传递
end

%% 3. 误比特率分析
for snrIndex = 1:numel(snrValues)
    totalErrorBits = 0;
    totalErrorPackets = 0;
    % 将Eb/No转换为SNR(考虑调制效率和码率)
    snrDb = snrValues(snrIndex) + 10*log10(log2(modulationOrder)*codeRate);

    for frameIndex = 1:frameCount
        % ===== 3.1 调制 =====
        switch MODULATION_TYPE
            case {0, 1, 2, 3}  % BPSK, QPSK, 16QAM, 64QAM
                % 比特到符号映射
                symbolBits = reshape(spreadBits, log2(modulationOrder), [])';  % 按调制阶数分组
                symbolIndices = bi2de(symbolBits, 'left-msb');  % 二进制转十进制索引
                
                % 基带调制
                if MODULATION_TYPE <= 1  % BPSK或QPSK
                    baseModulatedSignal = pskmod(symbolIndices, modulationOrder, pi/modulationOrder);
                else  % 16QAM或64QAM
                    baseModulatedSignal = qammod(symbolIndices, modulationOrder, 'UnitAveragePower', true);
                end
                
                % 上采样
                modulatedSignal = upsample(baseModulatedSignal, samplesPerBit);
                
                % 可选：脉冲成形滤波
                % span = 6;
                % rolloff = 0.25;
                % rrcFilter = rcosdesign(rolloff, span, samplesPerBit);
                % modulatedSignal = filter(rrcFilter, 1, modulatedSignal);
                
            case 4  % MSK
                % MSK调制(内部已包含上采样)
                bitsToModulate = spreadBits';
                modulatedSignal = mskModulator(bitsToModulate);
        end
        
        % 功率归一化
        modulatedSignal = modulatedSignal / sqrt(mean(abs(modulatedSignal).^2));
        
        % ===== 3.2 补零对齐 =====
        hopCount = ceil(length(modulatedSignal) / (samplesPerBit*bitsPerHop));  % 计算所需跳频点数
        paddingSampleCount = hopCount * samplesPerBit * bitsPerHop - length(modulatedSignal);
        modulatedSignal = [modulatedSignal; zeros(paddingSampleCount, 1)];  % 末尾补零对齐
        
        % 记录补零数量(用于解调时去除)
        exceedBitCount = paddingSampleCount;
        % 将信号重组为矩阵形式(每行对应一个跳频点)
        modulatedSignalMatrix = reshape(modulatedSignal, samplesPerBit*bitsPerHop, hopCount);
        modulatedSignalMatrix = modulatedSignalMatrix.';  % 转置: 行=跳频点, 列=采样点
        
        % ===== 3.3 跳频与信道传输 =====
        if ENABLE_FH == 0
            % 不跳频: 直接传输基带信号
            transmittedSignal = reshape(modulatedSignalMatrix.', numel(modulatedSignalMatrix), 1);
            receivedNoisy = awgn(transmittedSignal, snrDb);  % 添加高斯白噪声
            receivedNoisyMatrix = reshape(receivedNoisy, numel(receivedNoisy) / hopCount, hopCount);
            receivedBasebandMatrix = receivedNoisyMatrix.';
        else
            % 跳频: 生成随机跳频序列
            frequencyHoppingIndices = randi([1, frequencyPointCount], 1, hopCount);
            transmitFrequencyTable = carrierSequence(frequencyHoppingIndices);
            % 跳频调制
            frequencyHoppedSignalMatrix = FHmodulator(modulatedSignalMatrix, transmitFrequencyTable, samplingRate);

            % 信道传输
            transmittedSignal = reshape(frequencyHoppedSignalMatrix.', 1, numel(frequencyHoppedSignalMatrix));
            receivedNoisy = awgn(transmittedSignal, snrDb);  % 添加高斯白噪声

            % 重组为矩阵形式
            receivedNoisyMatrix = reshape(receivedNoisy, samplesPerBit*bitsPerHop, hopCount);
            receivedNoisyMatrix = receivedNoisyMatrix.';  % 行=跳频点, 列=采样点

            % 跳频解调
            receivedBasebandMatrix = FHdemodulator(receivedNoisyMatrix, transmitFrequencyTable, samplingRate);
        end

        % ===== 3.4 解调 =====
        receivedBasebandSignal = reshape(receivedBasebandMatrix.', 1, numel(receivedBasebandMatrix));

        % 根据调制方式区分处理
        if MODULATION_TYPE == 4  % MSK
            % MSK专用解调流程
            receivedSignal = reshape(receivedBasebandSignal, [], 1);
            
            % 使用MSK解调器
            allDemodulatedBits = mskDemodulator(receivedSignal);
            
            % 去除补零比特
            if exceedBitCount > 0
                demodulatedBitStream = allDemodulatedBits(1:end-exceedBitCount)';
            else
                demodulatedBitStream = allDemodulatedBits';
            end
        else
            % 其他调制方式解调流程
            % 去除补零样本
            if exceedBitCount > 0
                receivedData = receivedBasebandSignal(1:end-exceedBitCount);
            else
                receivedData = receivedBasebandSignal;
            end
            
            % 转为列向量
            receivedSignal = reshape(receivedData, [], 1);
            
            % 下采样
            downsampled = receivedSignal(1:samplesPerBit:end);
            
            % 解调
            if MODULATION_TYPE <= 1  % BPSK或QPSK
                demodulatedIndices = pskdemod(downsampled, modulationOrder, pi/modulationOrder);
            else  % 16QAM或64QAM
                demodulatedIndices = qamdemod(downsampled, modulationOrder, 'UnitAveragePower', true);
            end
            
            % 符号到比特转换
            demodulatedBits = de2bi(demodulatedIndices, log2(modulationOrder), 'left-msb');
            demodulatedBitStream = reshape(demodulatedBits', 1, []);
        end

        %% 解扩频
        if ENABLE_CCSK == 1
            ccsk_dec = CCSKDemodulator(demodulatedBitStream, ccskCode);  % CCSK解扩
        else
            ccsk_dec = double(demodulatedBitStream);  % 不使用CCSK时直接传递
        end

        %% 解交织
        if ENABLE_INTERLEAVING == 1
            deinterleavedBits = zeros(1, length(ccsk_dec));
            for k = 1:interleavingCount
                startIdx = ((k-1)*(interleavingRows*interleavingCols))+1;
                endIdx = k*(interleavingRows*interleavingCols);
                tempData = ccsk_dec(1, startIdx:endIdx);
                % 按行填充矩阵，按列读出实现解交织
                deinterleavedData = matdeintrlv(tempData, interleavingRows, interleavingCols);
                deinterleavedBits(1:k*interleavingRows*interleavingCols) = ...
                    horzcat(deinterleavedBits(1:(k-1)*interleavingRows*interleavingCols), deinterleavedData);
            end
        else
            deinterleavedBits = ccsk_dec;  % 不使用交织时直接传递
        end

        %% 信道译码
        switch CODING_TYPE
            case 0  % 不编码
                decodedBits = deinterleavedBits;
            case 1  % RS解码
                decodedBits = RSDecoder(deinterleavedBits, codewordLength, messageLength);
            case 2  % 卷积解码
                decodedBits = vitdec(deinterleavedBits, trellis, tracebackDepth, 'trunc', 'hard');  % Viterbi硬判决译码
        end

        %% CRC检错
        if ENABLE_CRC == 1
            [decodedData, crcError] = crcDetector(decodedBits');  % CRC检错(输入为列向量)
            crcErrorCount(snrIndex) = crcErrorCount(snrIndex) + sum(crcError);
            decodedData = decodedData';  % 转换回行向量
        else
            decodedData = decodedBits;  % 不使用CRC时直接传递
        end

        % 计算误比特数
        % 确保比较的两个向量长度相同
        minLength = min(length(inputBits), length(decodedData));
        errorBits = sum(inputBits(1:minLength) ~= decodedData(1:minLength));
        totalErrorBits = totalErrorBits + errorBits;

        % 计算误包率
        originalPacket = inputBits;
        decodedPacket = decodedData;
        originalPacketMatrix = reshape(originalPacket, messageBitCount/packetCount, packetCount);
        decodedPacketMatrix = reshape(decodedPacket, messageBitCount/packetCount, packetCount);
        packetErrors = sum(bsxfun(@ne, originalPacketMatrix, decodedPacketMatrix), 1) > 0;
        totalErrorPackets = totalErrorPackets + sum(packetErrors);
    end

    % 计算误比特率和误包率
    berValues(snrIndex) = totalErrorBits / (frameCount * totalBitCount);
    packetErrorRate(snrIndex) = totalErrorPackets / (frameCount * packetCount);
end

%% 误码性能分析
figure
semilogy(snrValues, berValues, '-*')
save MSK+RS+CCSK+FH berValues
hold on
grid on

xlabel('Eb/No (dB)')  % 横轴: 比特能量与噪声功率谱密度比
ylabel('误比特率')     % 纵轴: 误比特率(BER)
title(['Link16系统性能分析: ', figureLabel_Modulation, ' + ', figureLabel_Coding, ' + ', figureLabel_CCSK, ' + ', figureLabel_FH]);

toc
