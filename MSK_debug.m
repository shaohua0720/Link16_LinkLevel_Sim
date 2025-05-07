%% MSK调制解调调试脚本
clc; clear; close all;

%% 参数设置
samp = 20;  % 过采样倍数

%% 测试单比特
fprintf('测试单比特调制解调...\n');

% 测试比特0
bit0 = [0];
[~, sig0] = MSKmodulator(samp, bit0);
demod0 = MSKdemodulator(samp, sig0);
fprintf('比特0 -> 解调结果: %d (应为0)\n', demod0);

% 测试比特1
bit1 = [1];
[~, sig1] = MSKmodulator(samp, bit1);
demod1 = MSKdemodulator(samp, sig1);
fprintf('比特1 -> 解调结果: %d (应为1)\n', demod1);

%% 测试两个比特
fprintf('\n测试两个比特调制解调...\n');

% 测试00
bits00 = [0 0];
[~, sig00] = MSKmodulator(samp, bits00);
demod00 = MSKdemodulator(samp, sig00);
fprintf('比特00 -> 解调结果: %d%d (应为00)\n', demod00(1), demod00(2));

% 测试01
bits01 = [0 1];
[~, sig01] = MSKmodulator(samp, bits01);
demod01 = MSKdemodulator(samp, sig01);
fprintf('比特01 -> 解调结果: %d%d (应为01)\n', demod01(1), demod01(2));

% 测试10
bits10 = [1 0];
[~, sig10] = MSKmodulator(samp, bits10);
demod10 = MSKdemodulator(samp, sig10);
fprintf('比特10 -> 解调结果: %d%d (应为10)\n', demod10(1), demod10(2));

% 测试11
bits11 = [1 1];
[~, sig11] = MSKmodulator(samp, bits11);
demod11 = MSKdemodulator(samp, sig11);
fprintf('比特11 -> 解调结果: %d%d (应为11)\n', demod11(1), demod11(2));

%% 可视化调制信号
figure;
subplot(4,1,1);
t = (0:length(sig00)-1)/samp;
plot(t, real(sig00), 'b-', t, imag(sig00), 'r-');
title('MSK信号 - 比特00');
legend('I分量', 'Q分量');

subplot(4,1,2);
plot(t, real(sig01), 'b-', t, imag(sig01), 'r-');
title('MSK信号 - 比特01');
legend('I分量', 'Q分量');

subplot(4,1,3);
plot(t, real(sig10), 'b-', t, imag(sig10), 'r-');
title('MSK信号 - 比特10');
legend('I分量', 'Q分量');

subplot(4,1,4);
plot(t, real(sig11), 'b-', t, imag(sig11), 'r-');
title('MSK信号 - 比特11');
legend('I分量', 'Q分量');

%% 可视化解调过程
% 以比特10为例
figure;
t_bit = (0:samp-1)/samp;
t_full = (0:length(sig10)-1)/samp;

% 载波
subplot(3,2,1);
plot(t_full, cos(2*pi*0.25*t_full));
title('I载波');

subplot(3,2,2);
plot(t_full, sin(2*pi*0.25*t_full));
title('Q载波');

% 调制信号
subplot(3,2,3);
plot(t_full, real(sig10));
title('I分量');

subplot(3,2,4);
plot(t_full, imag(sig10));
title('Q分量');

% 解调信号
I_demod = real(sig10) .* cos(2*pi*0.25*t_full);
Q_demod = imag(sig10) .* sin(2*pi*0.25*t_full);

subplot(3,2,5);
plot(t_full, I_demod);
title('I解调');

subplot(3,2,6);
plot(t_full, Q_demod);
title('Q解调');

%% 详细分析偶数位解调
% 以比特01为例，分析第二位（偶数位）的解调
figure;
t_bit = (0:samp-1)/samp;
t_full = (0:length(sig01)-1)/samp;

% 提取第二个比特的Q分量
bit_idx = samp+1:2*samp;
Q_signal = imag(sig01);
Q_bit = Q_signal(bit_idx);

% 载波
subplot(3,1,1);
plot(t_bit, sin(2*pi*0.25*(samp+t_bit)));
title('Q载波 (第二个比特)');

% Q分量
subplot(3,1,2);
plot(t_bit, Q_bit);
title('Q分量 (第二个比特)');

% 相关波形
subplot(3,1,3);
Q_carrier = sin(2*pi*0.25*(samp+t_bit));
Q_corr = Q_bit .* Q_carrier;
plot(t_bit, Q_corr);
title('Q相关 (第二个比特)');
corr_value = sum(Q_corr .* sin(pi*t_bit/2));
fprintf('01的第二位Q相关值: %.4f (应为负值表示1)\n', corr_value);
