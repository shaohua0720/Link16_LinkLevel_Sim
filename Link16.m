%% Link16通信系统仿真
% 版本: v1.0.0
% 简易跳频系统误比特率统计分析
% 调制-->跳频-->信道-->解跳-->解调-->误码分析
tic
clc; clear; close all;

%% 系统参数设置
% 系统功能开关
ENABLE_FH = 0;                  % 0: 不跳频传输 1: 跳频传输
MODULATION_TYPE = 2;            % 0: BPSK, 1: QPSK, 2: 16QAM, 3: 64QAM, 4: MSK
ENABLE_CRC = 0;                 % 0: 不使用CRC编码, 1: 使用CRC编码
CODING_TYPE = 1;                % 0: 不编码, 1: RS(31,15)编码, 2: (2,1,7)卷积编码
ENABLE_INTERLEAVING = 1;        % 0: 不使用交织, 1: 使用交织
ENABLE_CCSK = 1;                % 0: 不使用CCSK扩频, 1: 使用CCSK扩频

% 初始化补零计数变量
exceedBitCount = 0;             % 初始化补零计数变量

% 基本系统参数
bitRate = 5e4;                  % 比特率: 50Kb/s
samplesPerBit = 20;             % 过采样倍数
samplingRate = samplesPerBit * bitRate;  % 采样率: 1MHz

% 跳频参数设置
hoppingRate = 1000;             % 跳频速率: 1000跳/秒
bitsPerHop = bitRate / hoppingRate;  % 每跳比特数目: 50比特/跳
frequencyBandwidth = 5e6;       % 跳频带宽: 5MHz
frequencyPointCount = floor(frequencyBandwidth / (bitRate*4));  % 跳频频点数目: 25个频点
frequencyInterval = frequencyBandwidth / frequencyPointCount;  % 频点间隔: 200kHz
frequencySequence = ((0:frequencyPointCount-1) - floor(frequencyPointCount/2)) * frequencyInterval;  % 跳频频点序列: 相对于中心频率的偏移
carrierFrequency = 3e6;         % 跳频中心频率: 3MHz
carrierSequence = carrierFrequency + frequencySequence;  % 发送时跳频频点序列: 绝对频率值

% 信道参数
dopplerShift = 300;             % 多普勒频偏: 300Hz
pathPowerProfile = [-1.0 -1.0 -1.0 0 0 0 -3.0 -5.0 -7.0];  % 多径功率谱: 单位dB
pathDelayProfile = [0 50 120 200 230 500 1600 2300 5000]*1e-8;  % 多径时延谱: 单位秒
rayleighChannel = comm.RayleighChannel('SampleRate', samplingRate, ...
    'PathDelays', pathDelayProfile, 'AveragePathGains', pathPowerProfile, ...
    'MaximumDopplerShift', dopplerShift, 'FadingTechnique', 'Sum of sinusoids');

% 设置显示标签
if ENABLE_FH == 1
    figureLabel_FH = '跳频传输';
else
    figureLabel_FH = '非跳频传输';
end

% 调制参数设置
switch MODULATION_TYPE
    case 0  % BPSK
        bitsPerSymbol = 1;
        figureLabel_Modulation = 'BPSK';
    case 1  % QPSK
        bitsPerSymbol = 2;
        figureLabel_Modulation = 'QPSK';
    case 2  % 16QAM
        bitsPerSymbol = 4;
        figureLabel_Modulation = '16QAM';
    case 3  % 64QAM
        bitsPerSymbol = 6;
        figureLabel_Modulation = '64QAM';
    case 4  % MSK
        bitsPerSymbol = 1;  % MSK是1比特/符号，与BPSK相同
        figureLabel_Modulation = 'MSK';
end
modulationOrder = 2^bitsPerSymbol;

% CRC参数设置
if ENABLE_CRC == 1
    crcBlockSize = 225;
    figureLabel_CRC = 'CRC(237,225)';  % 225位信息 + 12位CRC校验
else
    crcBlockSize = 1;  % 不使用CRC时，设为1以避免影响后续计算
    figureLabel_CRC = 'no CRC';
end

% 编码参数设置
rs_m = 5;  % RS码中每个符号的比特数: 5比特/符号
rs_n = 2^rs_m - 1;  % 编码后码字长度: 31符号
rs_t = 8;  % 能纠正的符号错误个数: 8符号
rs_k = rs_n - 2*rs_t;  % 信息段长度: 15符号

switch CODING_TYPE
    case 0  % 不编码
        figureLabel_Coding = '无纠错编码';
    case 1  % RS(31,15)编码
        figureLabel_Coding = 'RS(31,15)编码';
    case 2  % (2,1,7)卷积编码
        figureLabel_Coding = '(2,1,7)卷积编码';
end

% 信息参数设置
messageLength = rs_m * rs_k * bitsPerSymbol * crcBlockSize;  % 输入大小: 75比特 (RS码要求输入为符号整数倍)
syncBitCount = 0;                                            % 同步bit数目: 0比特
frameCount = 1;                                              % 传输帧数目: 1帧
packetCount = 3;                                             % 每一帧的包数目: 3包
bitsPerPacket = messageLength * bitsPerHop;                  % 每一包的比特数: 3750比特
messageBitCount = bitsPerPacket * packetCount;               % 有效消息bit数目: 11250比特
totalBitCount = syncBitCount + messageBitCount;              % 需要发送的bit总数: 11250比特

% 创建MSK调制器和解调器对象
if MODULATION_TYPE == 4  % MSK
    % 方案1: 完全使用MATLAB内置MSK调制解调器
    mskModulator = comm.MSKModulator('BitInput', true, ...
                                     'SamplesPerSymbol', samplesPerBit, ...
                                     'InitialPhaseOffset', 0);  % 明确设置初始相位偏移

    mskDemodulator = comm.MSKDemodulator('BitOutput', true, ...
                                         'SamplesPerSymbol', samplesPerBit, ...
                                         'InitialPhaseOffset', 0);  % 使用相同的初始相位偏移

    % 方案2: 完全使用自定义MSK调制解调器
    % 调制时使用自定义MSK调制器
    % [~, modulatedSignal] = MSKmodulator(samplesPerBit, reshape(paddedBits, hopCount, bitsPerHop));

    % 解调时使用对应的自定义MSK解调器
    % demodulatedBitStream = MSKdemodulator(samplesPerBit, receivedSignal);
end

%% 构造发送序列
syncBits = randi([0, 1], 1, syncBitCount);                   % 同步二进制序列（用一串随机序列代替）
messageBits = randi([0, 1], 1, totalBitCount - syncBitCount);  % 随机生成消息比特
inputBits = [syncBits, messageBits];                          % 构造整个发送bit序列

%% CRC检错编码
if ENABLE_CRC
    polynomial = 'z12+1';     % CRC多项式: z^12+1
    checksumPerFrame = length(inputBits)/crcBlockSize;  % 将传入帧细分为checksumPerFrame个等长的子帧
    crcGenerator = comm.CRCGenerator(polynomial, 'ChecksumsPerFrame', checksumPerFrame);  % CRC编码生成器
    crcDetector = comm.CRCDetector(polynomial, 'ChecksumsPerFrame', checksumPerFrame);  % CRC解码检测器
    crcEncodedBits = crcGenerator(inputBits');  % CRC编码 输入应是列向量
else
    crcEncodedBits = inputBits';  % 不使用CRC时直接传递
end
crcEncodedBits = crcEncodedBits';  % 转换回行向量

%% 信源纠错编码
switch CODING_TYPE
    case 0  % 不编码
        codedBits = crcEncodedBits;
    case 1  % RS编码
        % (31,15)RS编码（每15个符号用31个符号来表示）
        codewordLength = rs_n;  % 编码后码字长度: 31符号
        messageLength = rs_k;   % 信息段长度: 15符号
        codedBits = RSEncoder(crcEncodedBits, codewordLength, messageLength);  % RS编码
    case 2  % 卷积编码
        % (2,1,7)卷积码: 码率1/2，约束长度7
        constraintLength = 7;  % 约束长度: 7
        tracebackDepth = 42;   % Viterbi译码器回溯深度: 42 (约束长度的6倍)
        trellis = poly2trellis(constraintLength, [171, 133]);  % 卷积码生成多项式: 171和133 (八进制)
        codedBits = convenc(crcEncodedBits, trellis);  % 卷积编码
end

%% 交织 (对原数据分块进行处理：只是改变数据位置 不改变数据数量)
if ENABLE_INTERLEAVING
    interleavingRows = 10;
    interleavingCols = 100;  % 设定交织的深度与宽度: 10×100矩阵
    interleavingCount = length(codedBits)/(interleavingRows*interleavingCols);  % 交织次数
    interleavedBits = zeros(1, length(codedBits));

    for i = 1:interleavingCount
        startIdx = ((i-1)*(interleavingRows*interleavingCols))+1;
        endIdx = i*(interleavingRows*interleavingCols);
        tempData = codedBits(1, startIdx:endIdx);
        % 通过按列填充矩阵，并按行输出来实现交织
        interleavedData = matintrlv(tempData, interleavingRows, interleavingCols);
        interleavedBits(1:i*interleavingRows*interleavingCols) = ...
            horzcat(interleavedBits(1:(i-1)*interleavingRows*interleavingCols), interleavedData);
    end
    figureLabel_Interleaving = '使用交织';
else
    interleavedBits = codedBits;
    figureLabel_Interleaving = '不使用交织';
end

%% CCSK软扩频 CCSK序列为32位
if ENABLE_CCSK == 1
    % 优化后的CCSK序列 相关性更好
    ccskCode = [1 0 1 1 1 0 1 0 0 0 1 1 1 1 0 1 0 0 1 0 0 0 0 0 0 1 1 0 0 1 1 0]';
    spreadBits = CCSKModulator(interleavedBits, ccskCode);  % CCSK扩频: 每1比特扩展为32比特
    figureLabel_CCSK = 'CCSK';
else
    spreadBits = interleavedBits;
    figureLabel_CCSK = '不使用CCSK';
end

%% 星座映射
symbolBits = reshape(spreadBits, log2(modulationOrder), [])';  % 以每组log2(M)比特进行分组
symbolIndices = bi2de(symbolBits, 'left-msb');  % 二进制转化为十进制

switch MODULATION_TYPE
    case 0  % BPSK
        modulatedSignal = pskmod(symbolIndices, modulationOrder, pi/modulationOrder);  % 已经归一化
    case 1  % QPSK
        modulatedSignal = pskmod(symbolIndices, modulationOrder, pi/modulationOrder);  % 已经归一化
    case 2  % 16QAM
        modulatedSignal = qammod(symbolIndices, modulationOrder, 'UnitAveragePower', true);
    case 3  % 64QAM
        modulatedSignal = qammod(symbolIndices, modulationOrder, 'UnitAveragePower', true);
    case 4  % MSK
        % MSK专用处理流程
        % 1. 不分割为跳频点，而是整体调制
        bitsToModulate = spreadBits';

        % 2. 调制整个比特流
        modulatedSignal = mskModulator(bitsToModulate);

        % 3. 功率归一化
        modulatedSignal = modulatedSignal / sqrt(mean(abs(modulatedSignal).^2));

        % 4. 然后再分割为跳频点
        hopCount = ceil(length(modulatedSignal) / (samplesPerBit*bitsPerHop));
        totalSamplesNeeded = hopCount * samplesPerBit * bitsPerHop;

        % 5. 样本级别补零
        paddingSampleCount = 0;  % 初始化补零计数变量
        if length(modulatedSignal) < totalSamplesNeeded
            paddingSampleCount = totalSamplesNeeded - length(modulatedSignal);
            modulatedSignal = [modulatedSignal; zeros(paddingSampleCount, 1)];
        end

        % 6. 重塑为矩阵形式
        modulatedSignalMatrix = reshape(modulatedSignal, samplesPerBit*bitsPerHop, hopCount).';

        % 7. 记录补零信息
        exceedBitCount = paddingSampleCount;
end

averagePower = mean(abs(modulatedSignal).^2);  % 计算信号平均功率
% scatterplot(modulatedSignal); title([figureLabel_Modulation, '星座图']);

%% 补零(使数据长度为跳频点数的整数倍)
hopCount = ceil(length(modulatedSignal) / (samplesPerBit*bitsPerHop));  % 发送所有的bit需要的跳频点数
% 由于总的发送bit数目不一定是每跳bit数目bitsPerHop的整数倍，所以，构造为bitsPerHop的整数倍，在解调后去除多余bit
paddingSampleCount = hopCount * samplesPerBit * bitsPerHop - length(modulatedSignal);
modulatedSignal = [modulatedSignal; zeros(paddingSampleCount, 1)];  % 在信号末尾补零
modulatedSignalMatrix = reshape(modulatedSignal, samplesPerBit*bitsPerHop, hopCount);  % 将待发送序列转成一个矩阵
modulatedSignalMatrix = modulatedSignalMatrix.';  % 转置: 行对应跳频点，列对应每跳的采样点

%% 误比特分析参数设置
snrValues = (-10:0);  % Eb/No值范围: -10dB到0dB，步长1dB
codeRate = 1/2;  % 码率: 1/2 (考虑了编码和扩频的影响)
berValues = zeros(1, length(snrValues));  % 误比特率数组
packetErrorRate = zeros(1, length(snrValues));  % 误包率数组
crcErrorCount = zeros(1, length(snrValues));  % CRC检测到的错误计数

%% 误比特率分析
for snrIndex = 1:numel(snrValues)
    totalErrorBits = 0;
    totalErrorPackets = 0;
    snrDb = snrValues(snrIndex) + 10*log10(log2(modulationOrder)*codeRate);  % 将Eb/No转换为SNR

    for frameIndex = 1:frameCount
        %% 跳频
        if ENABLE_FH == 0
            transmittedSignal = reshape(modulatedSignalMatrix.', numel(modulatedSignalMatrix), 1);  % 将矩阵形式的信号转化为实际的1维信号
            receivedNoisy = awgn(transmittedSignal, snrDb);  % 添加高斯白噪声
            receivedNoisyMatrix = reshape(receivedNoisy, numel(receivedNoisy) / hopCount, hopCount);  % 为了便于后面的解跳/解跳操作, 将1维信号还原为矩阵形式
            receivedBasebandMatrix = receivedNoisyMatrix.';
        else
            frequencyHoppingIndices = randi([1, frequencyPointCount], 1, hopCount);  % 根据跳频点数生成随机频点序列索引
            transmitFrequencyTable = carrierSequence(frequencyHoppingIndices);  % 根据随机频点序列索引生成跳频频点（收发频点保持一致才可以解跳）
            frequencyHoppedSignalMatrix = FHmodulator(modulatedSignalMatrix, transmitFrequencyTable, samplingRate);  % 跳频调制

            %% 信道
            transmittedSignal = reshape(frequencyHoppedSignalMatrix.', 1, numel(frequencyHoppedSignalMatrix));  % 将矩阵形式的信号转化为实际的1维信号
            receivedNoisy = awgn(transmittedSignal, snrDb);  % 添加高斯白噪声

            receivedNoisyMatrix = reshape(receivedNoisy, samplesPerBit*bitsPerHop, hopCount);  % 使用与发送时相同的维度
            receivedNoisyMatrix = receivedNoisyMatrix.';  % 行：跳数，列: 每跳对应的带噪声的调制信号

            %% 解跳并去掉冗余bit
            receivedBasebandMatrix = FHdemodulator(receivedNoisyMatrix, transmitFrequencyTable, samplingRate);  % 跳频解调
        end

        receivedBasebandSignal = reshape(receivedBasebandMatrix.', 1, numel(receivedBasebandMatrix));

        % 根据调制方式区分处理
        if MODULATION_TYPE == 4  % MSK
            % MSK专用处理流程
            % 1. 保留完整信号用于MSK解调
            receivedSignal = reshape(receivedBasebandSignal, [], 1);

            % 2. 相位校正（可选）
            % 估计相位偏移
            phaseOffset = angle(mean(receivedSignal));
            % 校正相位
            receivedSignalCorrected = receivedSignal * exp(-1j * phaseOffset);

            % 3. 使用MATLAB内置MSK解调器
            allDemodulatedBits = mskDemodulator(receivedSignalCorrected);

            % 4. 去除补零比特
            if exceedBitCount > 0
                demodulatedBitStream = allDemodulatedBits(1:end-exceedBitCount)';
            else
                demodulatedBitStream = allDemodulatedBits';
            end
        else
            % 其他调制方式处理流程
            % 1. 样本级别去零
            if exceedBitCount > 0
                receivedData = receivedBasebandSignal(1:end-exceedBitCount);
            else
                receivedData = receivedBasebandSignal;
            end

            % 2. 重塑为列向量
            receivedSignal = reshape(receivedData, [], 1);

            % 3. 根据调制方式解调
            switch MODULATION_TYPE
                case 0  % BPSK
                    demodulatedIndices = pskdemod(receivedSignal, modulationOrder, pi/modulationOrder);
                case 1  % QPSK
                    demodulatedIndices = pskdemod(receivedSignal, modulationOrder, pi/modulationOrder);
                case 2  % 16QAM
                    demodulatedIndices = qamdemod(receivedSignal, modulationOrder, 'UnitAveragePower', true);
                case 3  % 64QAM
                    demodulatedIndices = qamdemod(receivedSignal, modulationOrder, 'UnitAveragePower', true);
            end

            % 4. 转换为比特流
            demodulatedIndices = reshape(demodulatedIndices, [], 1);
            demodulatedBits = de2bi(demodulatedIndices, log2(modulationOrder), 'left-msb');
            demodulatedBitStream = reshape(demodulatedBits', 1, []);
        end

        %% 解扩
        if ENABLE_CCSK == 1
            ccsk_dec = CCSKDemodulator(demodulatedBitStream, ccskCode);
        else
            ccsk_dec = double(demodulatedBitStream);
        end

        %% 解交织
        if ENABLE_INTERLEAVING == 1
            deinterleavedBits = zeros(1, length(ccsk_dec));
            for k = 1:interleavingCount
                startIdx = ((k-1)*(interleavingRows*interleavingCols))+1;
                endIdx = k*(interleavingRows*interleavingCols);
                tempData = ccsk_dec(1, startIdx:endIdx);
                deinterleavedData = matdeintrlv(tempData, interleavingRows, interleavingCols);
                deinterleavedBits(1:k*interleavingRows*interleavingCols) = ...
                    horzcat(deinterleavedBits(1:(k-1)*interleavingRows*interleavingCols), deinterleavedData);
            end
        else
            deinterleavedBits = ccsk_dec;
        end

        %% 信道译码
        switch CODING_TYPE
            case 0  % 'NO'
                decodedBits = deinterleavedBits;
            case 1  % 'rs'
                decodedBits = RSDecoder(deinterleavedBits, codewordLength, messageLength);
            case 2  % 'juanji'
                decodedBits = vitdec(deinterleavedBits, trellis, tracebackDepth, 'trunc', 'hard');  % Viterbi译码
        end

        %% 检错解码
        if ENABLE_CRC == 1
            [decodedData, crcError] = crcDetector(decodedBits');  % 输入为列向量
            crcErrorCount(snrIndex) = crcErrorCount(snrIndex) + sum(crcError);
            decodedData = decodedData';
        else
            decodedData = decodedBits;
        end

        % 计算误比特数
        % 确保比较的两个向量长度相同
        minLength = min(length(inputBits), length(decodedData));
        errorBits = sum(inputBits(1:minLength) ~= decodedData(1:minLength));
        totalErrorBits = totalErrorBits + errorBits;

        % 计算误包率
        originalPacket = inputBits(syncBitCount+1:end);
        decodedPacket = decodedData(syncBitCount+1:end);
        originalPacketMatrix = reshape(originalPacket, messageBitCount/packetCount, packetCount);
        decodedPacketMatrix = reshape(decodedPacket, messageBitCount/packetCount, packetCount);
        packetErrors = sum(bsxfun(@ne, originalPacketMatrix, decodedPacketMatrix), 1) > 0;
        totalErrorPackets = totalErrorPackets + sum(packetErrors);
    end

    % 计算误比特率和误包率
    berValues(snrIndex) = totalErrorBits / (frameCount * totalBitCount);
    packetErrorRate(snrIndex) = totalErrorPackets / (frameCount * packetCount);
end

%% 误码输出
figure
semilogy(snrValues, berValues, '-*')
save MSK+RS+CCSK+FH berValues
hold on
grid on

xlabel('Eb/No (dB)')  % 如果不进行转换就是SNR
ylabel('误比特率')     % 如果不进行转换就是误码率
title(['Link16系统性能分析: ', figureLabel_Modulation, ' + ', figureLabel_Coding, ' + ', figureLabel_CCSK, ' + ', figureLabel_FH]);
legend('MSK+RS+CCSK+FH')

toc
