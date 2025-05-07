%% MSK调制解调测试脚本
clc; clear; close all;

%% 参数设置
samp = 20;                % 过采样倍数
bitCount = 1000;          % 测试比特数
snrValues = -10:2:10;     % 测试的信噪比范围

%% 生成随机测试数据
txBits = randi([0, 1], 1, bitCount);

%% 单跳频点测试
TX_BIT_MAT = txBits;      % 单跳频点情况下直接使用比特序列

%% MSK调制
[biNRZ, txSig] = MSKmodulator(samp, TX_BIT_MAT);

%% 在不同信噪比下测试
berValues = zeros(1, length(snrValues));

for i = 1:length(snrValues)
    snr = snrValues(i);
    
    % 添加噪声
    rxSig = awgn(txSig, snr);
    
    % MSK解调
    demodBits = MSKdemodulator(samp, rxSig);
    
    % 计算误比特率
    errors = sum(txBits ~= demodBits(1:length(txBits)));
    berValues(i) = errors / bitCount;
    
    % 显示进度
    fprintf('SNR = %d dB, BER = %.6f\n', snr, berValues(i));
end

%% 绘制误比特率曲线
figure;
semilogy(snrValues, berValues, '-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('误比特率 (BER)');
title('MSK调制解调性能测试');

%% 理论MSK误比特率曲线
theoreticalBer = 0.5 * erfc(sqrt(10.^(snrValues/10)/2));
hold on;
semilogy(snrValues, theoreticalBer, '--r', 'LineWidth', 2);
legend('实际性能', 'MSK理论性能');

%% 波形可视化测试
% 选择一小段数据进行可视化
visualBits = txBits(1:20);
[~, visualTxSig] = MSKmodulator(samp, visualBits);

% 解调
visualDemodBits = MSKdemodulator(samp, visualTxSig);

% 绘制原始比特、调制信号和解调比特
figure;
subplot(3,1,1);
stem(visualBits, 'filled');
title('原始比特序列');
ylim([-0.2, 1.2]);
grid on;

subplot(3,1,2);
t = (0:length(visualTxSig)-1)/samp;
plot(t, real(visualTxSig), 'b-', t, imag(visualTxSig), 'r-');
title('MSK调制信号');
legend('I分量', 'Q分量');
grid on;

subplot(3,1,3);
stem(visualDemodBits(1:length(visualBits)), 'filled');
title('解调比特序列');
ylim([-0.2, 1.2]);
grid on;

%% 多跳频点测试
hopCount = 5;
bitsPerHop = 200;
TX_BIT_MAT_MULTI = randi([0, 1], hopCount, bitsPerHop);

% MSK调制
[biNRZ_multi, txSig_multi] = MSKmodulator(samp, TX_BIT_MAT_MULTI);

% 无噪声情况下解调
for hop = 1:hopCount
    rxSig_hop = txSig_multi(hop, :);
    demodBits_hop = MSKdemodulator(samp, rxSig_hop);
    
    % 计算每个跳频点的误比特率
    errors = sum(TX_BIT_MAT_MULTI(hop, :) ~= demodBits_hop(1:bitsPerHop));
    fprintf('跳频点 %d: 误比特率 = %.6f\n', hop, errors/bitsPerHop);
end

%% 相位连续性测试
% 生成连续的比特序列
contBits = randi([0, 1], 1, 50);
[~, contTxSig] = MSKmodulator(samp, contBits);

% 绘制相位连续性
figure;
t = (0:length(contTxSig)-1)/samp;
plot(t, angle(contTxSig));
title('MSK信号相位连续性');
xlabel('时间 (比特周期)');
ylabel('相位 (弧度)');
grid on;