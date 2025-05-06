%% 知乎链接:https://www.zhihu.com/people/san-hao-bai-du-ren-79
%% 微信公众号：一个安静的资料号

%% 简易跳频系统误比特率统计分析
% 调制-->跳频-->信道-->解跳-->解调-->误码分析
tic
clc;clear;close all;

%% 参数设置
FH = 1;                         % 0: 不跳频传输 1：跳频传输
XINGZUO = 4;%BPSK:0   QPSK:1  '16QAM':2   '64QAM':3    'MSK':4
CRC_CODE = 0;%1： CRC编码  0：不CRC编码
DatainCode = 1;%'NO': 0 不编码 || 'rs'：1 RS(31,15)编码  || 'juanji'：2 (2,1,7)卷积编码 
interwine = 1;
CCSK = 1;
% ******************** 系统参数设置*****************%
Rb = 5e4;                       % 速率：50Kb/s
Tb = 1 / Rb;                    % bit间隔

% ********************跳频参数设置 *****************%
hopping = 1000;                 % 跳频速率
bitsPerHop = Rb / hopping;      % 每跳bit数目（必须为整数）
samp = 20;                      % 过采样倍数
fs = samp*Rb;                   % 采样率
BW = 5e6;                       % 跳频带宽
freqNum = floor(BW / (Rb*4));   % 跳频频点数目
freqInterval = BW / freqNum;    % 频点间隔
freqSeq = ([0:freqNum-1] - floor(freqNum/2)) * freqInterval;   % 跳频频点序列
carrier = 3e6;                  % 跳频中心频率
carrierSeq = carrier + freqSeq; % 发送时跳频频点序列

ts = 1/fs;
fd = 300;                                          %多普勒频偏
pathPower = [-1.0 -1.0 -1.0 0 0 0 -3.0 -5.0 -7.0];
pathDelays = [0 50 120 200 230 500 1600 2300 5000]*1e-8;
rchan = comm.RayleighChannel('SampleRate',fs, ...
    'PathDelays',pathDelays,'AveragePathGains',pathPower, ...
    'MaximumDopplerShift',fd,'FadingTechnique','Sum of sinusoids');

switch FH
    case 0
        FIG_FH = '非跳频传输';
    case 1
        FIG_FH = '跳频传输';
end
% ********************调制参数设置 *****************%
switch XINGZUO
    case 0  %'BPSK'
        xz_n = 1;
        FIG_XINGZUO = 'BPSK';
    case 1  %'QPSK'
        xz_n = 2;
        FIG_XINGZUO = 'QPSK';
    case 2  %'16QAM'
        xz_n = 4;
        FIG_XINGZUO = '16QAM';
    case 3  %'64QAM'
        xz_n = 6;
        FIG_XINGZUO = '64QAM';
    case 4  %'MSK'
        xz_n = 4;
        FIG_XINGZUO = 'MSK';
end
M = 2^xz_n;

% ***************CRC(237,225) X^12+1 *****************%
switch CRC_CODE
    case 0  %不采取检错编码
        crc_n = 1;
    case 1  %(237,225)检错编码
        crc_n = 225;
end
% ***************RS(31,15)编码参数*****************%
switch DatainCode
    case 0  %0：不编码
%         rs_m = 1;rs_k = 1;
        rs_m = 5; rs_n = 2^rs_m - 1;
        rs_t = 8; rs_k = rs_n - 2*rs_t ;%t: 能纠正的符号错误个数 k: 信息段长度
        FIG_code = '无纠错编码';
    case 1  %1： RS(31,15)编码
        rs_m = 5; rs_n = 2^rs_m - 1;
        rs_t = 8; rs_k = rs_n - 2*rs_t ;%t: 能纠正的符号错误个数 k: 信息段长度
        FIG_code = 'RS(31,15)编码';
    case 2  %2: (2,1,7)卷积编码 
%         rs_m = 1;rs_k = 1;
        rs_m = 5; rs_n = 2^rs_m - 1;
        rs_t = 8; rs_k = rs_n - 2*rs_t ;%t: 能纠正的符号错误个数 k: 信息段长度
        FIG_code = '(2,1,7)卷积编码';
end
% ************************产生信源信号*******************%
MSG_len = rs_m * rs_k * xz_n*crc_n; % 输入大小最好为rs_m * rs_k*xz_n的整数倍(后面会每rs_m位转换为十进制，然后gf的输入的列数必须为rs_k)
% ***************** 传输信息参数设置*****************%
SYNC_BIT_NUM = 0;                                            % 同步bit数目
frameNum = 1;                                                % 传输帧数目
PACKET_NUM = 3;                                                % 每一帧的包数目
BIT_PER_PACKET = MSG_len*bitsPerHop;                                       %每一包的比特数
MSG_BIT_NUM = BIT_PER_PACKET*PACKET_NUM;                      % 有效消息bit数目
TX_BIT_NUM = SYNC_BIT_NUM + MSG_BIT_NUM;            % 需要发送的bit数目

%% 构造发送序列
SYNC = randi([0 ,1] , 1 , SYNC_BIT_NUM);                        % 同步二进制序列（用一串随机序列代替）
MSG = randi([0,1] ,1 ,TX_BIT_NUM - SYNC_BIT_NUM);       % 消息字符号
dataIn = [SYNC , MSG];                                          % 构造整个发送bit序列                                      % 

%% CRC检错编码
switch CRC_CODE
    case 0
        dataIn_crc = dataIn';
        FIG_CRC_CODE = 'no CRC';
    case 1  %(237,225)检错编码
        poly = 'z12+1';     %多项式
        ChecksumsPerFrame = length(dataIn)/crc_n;%将传入帧细分为ChecksumsPerFrame个等长的子帧。因为每crc_n位数据编码一次
        crcgenerator = comm.CRCGenerator(poly,'ChecksumsPerFrame',ChecksumsPerFrame);%CRC编码生成器
        crcdetector = comm.CRCDetector(poly,'ChecksumsPerFrame',ChecksumsPerFrame);%CRC解码生成器

        dataIn_crc = crcgenerator(dataIn'); %CRC编码 输入应是列向量
        FIG_CRC_CODE = 'CRC(237,225)';
end
dataIn_crc = dataIn_crc';%转换回行向量

%% 信源纠错编码
switch DatainCode
    case 0
        encoded_msg = dataIn_crc;                      %1*30 0000
    case 1  %'rs'  每kk位编码为nn位
        %下面是(31,15)RS编码（每15个用31个来表示）
        nn = rs_n;                        %编码后码字长度(信息段+监督段)
        kk = rs_k;                        %信息段长度
        encoded_msg = LSY_RSCode(dataIn_crc,nn,kk);%1 * 62 0000
        %encoded_msg2 = LSY_RSCode(dataIn_crc,nn,kk);%第二种RS编码写法 这两种都是对的
    case 2  %'juanji' 每k位编码为n位
        %(n,k,m) n为输出长度 k为输入长度 m为编码约束度 n*(m+1)为约束长度 k/n为码率
        %下面是(2,1,7)卷积码
        L = 7;                   %约束长度
        tbdepth = 42;       % Traceback depth for Viterbi decoder Viterbi译码器回溯深度 一般是约束长度的5-9倍
        trel = poly2trellis(L,[171,133]);         %卷积码生成多项式
        encoded_msg = convenc(dataIn_crc,trel);       %卷积编码 1*60 0000
end

%% 交织 (对原数据分块进行处理：只是改变数据位置 不改变数据数量)
interwine_msg = zeros(1,length(encoded_msg));
switch interwine
    case 0
        interwine_msg = encoded_msg;
        FIG_interwine = 'NO interwine';
    case 1
        rows = 10;cols = 100;%设定交织的深度与宽度 rows*cols应该等于输入矩阵的行数(因为交织按列取来做填充)
        division = length(encoded_msg)/(rows*cols);%交织次数
        for i =1:division
            temp_data_1 = encoded_msg(1,(((i-1)*(rows*cols))+1):(i*(rows*cols)));%通过按列填充矩阵，并按行输出符号来恢复符号排序 https://www.jianshu.com/p/ac6c18fc3545
            temp_data_2 = matintrlv(temp_data_1,rows,cols);
            interwine_msg(1:i*rows*cols) = horzcat(interwine_msg(1:(i-1)*rows*cols),temp_data_2);
        end
        FIG_interwine = 'interwine';
end
%% CCSK软扩频 CCSK序列为32位

switch CCSK
    case 0
        ccsk_msg = interwine_msg;
        FIG_CCSK = 'NO CCSK';
    case 1
        % ccskcode=[0 1 1 1 1 1 0 0 1 1 1 0 1 0 0 1 0 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0]';%Link16的CCSK序列 自相关性并不是最优
        ccskcode=[1 0 1 1 1 0 1 0 0 0 1 1 1 1 0 1 0 0 1 0 0 0 0 0 0 1 1 0 0 1 1 0]';%优化后的CCSK序列 相关性更好  Link16 关键技术研究及系统建模
        %n_ccsk = 32;%CCSK序列长度
        ccsk_msg = LSY_CCSK32(interwine_msg,ccskcode);
        FIG_CCSK = 'NO CCSK';
end
%% 星座映射
data_xz= reshape(ccsk_msg,log2(M),[])';       %以每组log2(M)比特进行分组，M=4
data_xzde= bi2de(data_xz,'left-msb');            %二进制转化为十进制

switch XINGZUO
    case 0  %'BPSK'
        txSig = pskmod(data_xzde,M,pi/M);               %已经归一化
    case 1  %'QPSK'
        txSig = pskmod(data_xzde,M,pi/M);               %已经归一化
    case 2  %'16QAM'
        txSig = qammod(data_xzde,M,'UnitAveragePower',true);
    case 3  %'64QAM'
        txSig = qammod(data_xzde,M,'UnitAveragePower',true);
    case 4  %'MSK'
        TX_BIT_MAT = reshape(ccsk_msg , bitsPerHop , []);                        %如果是MSK 则按频点处理
        TX_BIT_MAT = TX_BIT_MAT';  
        [biNRZ , txSig0]  = MSKmodulator(samp, TX_BIT_MAT); % 复基带调制矩阵
        txSig = reshape(txSig0' , [] , 1);
end

avgPower = mean(abs(txSig).^2);
% scatterplot(txSig);title([FIG_XINGZUO,'星座图 ']);
    
%% 补零(因为经过一系列处理后 数据不再是整数矩阵)
HOP_NUM = ceil(length(txSig) / (samp*bitsPerHop));            % 发送所有的bit需要的跳频点数
%（由于总的发送bit数目不一定是每跳bit数目bitsPerHop的整数倍，所以，构造为bitsPerHop的整数倍，在解调后去除多余bit即可）
EXCEED_BIT = HOP_NUM * samp * bitsPerHop - length(txSig);     % 多余的bit 
txSig = [txSig;zeros(EXCEED_BIT,1)];
% msgModMatrix = reshape(txSig , HOP_NUM , samp*bitsPerHop);            % 将待发送序列转成一个矩阵，行：跳频点数，列: 每跳对应的bit序列
msgModMatrix = reshape(txSig , samp*bitsPerHop , HOP_NUM);            % 将待发送序列转成一个矩阵，行：跳频点数，列: 每跳对应的bit序列
msgModMatrix = msgModMatrix.';

%% 误比特分析参数设置
 snr = (-10:1:10);                         % Eb/No values (dB)
%snr = (0:0.2:15);                         % Eb/No values (dB)
rate = 1/2;                             % 码率
ber = zeros(1 , length(snr));           % 误bit率
berPacket = zeros(1 , length(snr));     % 误包率
CRC = zeros(1,length(snr));

%% 误bit率分析（通过在不同的信噪比下传输多个帧做统计分析，可以通过增加传输帧数目获取更加精确的误比特率；更改信噪比获取不同信噪比下的性能）
for ii = 1 : length(snr)
    sumErrBit = 0;
    sumErrPacket = 0;
    snrdB = snr(ii) + 10*log10(log2(M)*rate);   % Convert Eb/No to SNR 
    for jj = 1 : frameNum
        %% 跳频
switch FH
    case 0
        txFHmodulated = reshape(msgModMatrix.' , numel(msgModMatrix) , 1);       % 将矩阵形式的信号转化为实际的1维信号 
        N_Interval = 200;N_Length= 4;
%         txFHmodulated = tufaxindao(N_Interval,N_Length,txFHmodulated);
%         channel_out = step(rchan,txFHmodulated);
        rcvNoisy = awgn(txFHmodulated , snrdB);     % 添加噪声
%       rcvNoisyMat = reshape(txFHmodulated , numel(txFHmodulated) / HOP_NUM , HOP_NUM);  % 为了便于后面的解跳/解跳操作, 将1维信号还原为矩阵形式
        rcvNoisyMat = reshape(rcvNoisy , numel(rcvNoisy) / HOP_NUM , HOP_NUM);  % 为了便于后面的解跳/解跳操作, 将1维信号还原为矩阵形式
        rcvBBmat = rcvNoisyMat.';                 
    case 1
        fhIndex = randi([1 , freqNum] , 1 , HOP_NUM);                           % 根据跳频点数生成随机频点序列索引
        txFHtable = carrierSeq(fhIndex);                                        % 根据随机频点序列索引生成跳频频点（收发频点保持一致才可以解跳）
        txFHmodulatedMat = FHmodulator(msgModMatrix , txFHtable , fs);	% 跳频后的信号
        
        %% 信道
        txFHmodulated = reshape(txFHmodulatedMat.' , 1 , numel(txFHmodulatedMat));       % 将矩阵形式的信号转化为实际的1维信号
        rcvNoisy = awgn(txFHmodulated , snrdB);     % 添加噪声
        
        rcvNoisyMat = reshape(rcvNoisy , numel(rcvNoisy) / HOP_NUM , HOP_NUM);  % 为了便于后面的解跳/解跳操作, 将1维信号还原为矩阵形式
        rcvNoisyMat = rcvNoisyMat.';                                                     % 行：跳数，列: 每跳对应的带噪声的调制信号
        % scatterplot(rxSig)
        % eyediagram(z,16);%眼图
        %% 解跳并去掉冗余bit
        rcvBBmat = FHdemodulator(rcvNoisyMat , txFHtable , fs);                                 % 根据冗余bit数目去除多余bit
end
        rcvBB1dim = reshape(rcvBBmat.' , 1, numel(rcvBBmat));
        rcv_data = rcvBB1dim(1:end-EXCEED_BIT);    %去掉末尾的0 

        %% 星座解调               % 还原信号为连续形式
        rxSig = reshape(rcv_data,[],1);
        switch XINGZUO
            case 0  %'BPSK'
                data_dec = pskdemod(rxSig,M,pi/M);
                data_dec = reshape(data_dec,[],1);
                data_dec_bi = de2bi(data_dec, log2(M), 'left-msb');%转换成m位二进制
                dmdBit = reshape(data_dec_bi',1,[]);
            case 1  %'QPSK'
                data_dec = pskdemod(rxSig,M,pi/M);
                data_dec = reshape(data_dec,[],1);
                data_dec_bi = de2bi(data_dec, log2(M), 'left-msb');%转换成m位二进制
                dmdBit = reshape(data_dec_bi',1,[]);
            case 2  %'16QAM'
                data_dec = qamdemod(rxSig,M,'UnitAveragePower',true);
                data_dec = reshape(data_dec,[],1);
                data_dec_bi = de2bi(data_dec, log2(M), 'left-msb');%转换成m位二进制
                dmdBit = reshape(data_dec_bi',1,[]);
            case 3  %'64QAM'
                data_dec = qamdemod(rxSig,M,'UnitAveragePower',true);
                data_dec = reshape(data_dec,[],1);
                data_dec_bi = de2bi(data_dec, log2(M), 'left-msb');%转换成m位二进制
                dmdBit = reshape(data_dec_bi',1,[]);
            case 4  %'MSK'
                % 差分解调
                rxSig = reshape(rxSig , 1 , numel(rxSig));                   % 还原信号为连续形式
                rcvBBsamp = reshape(rxSig , samp , numel(rxSig)/samp);          % 转换成每samp一个单元(即以bit为采样单元)
                dmdBit = imag(conj(rcvBBsamp(samp,:)) .* rcvBBsamp(1,:)) >0;            %conj函数:用于计算复数的共轭值 
        end
        
%% 解扩
switch CCSK
    case 0
        ccsk_dec = double(dmdBit);
    case 1
        ccsk_dec = LSY_CCSKde32(dmdBit,ccskcode);
end
%% 解交织
deinterwine_data = zeros(1,length(ccsk_dec));
switch interwine
    case 0
        deinterwine_data = ccsk_dec;
    case 1
        for k = 1:division
            temp_data_1 = ccsk_dec(1,(((k-1)*(rows*cols))+1):(k*(rows*cols)));
            temp_data_2 = matdeintrlv(temp_data_1,rows,cols);
            deinterwine_data(1:k*rows*cols) = horzcat(deinterwine_data(1:(k-1)*rows*cols),temp_data_2);
        end
end
%% 信道译码
switch DatainCode
    case 0  %'NO'
        dataOut_dc = deinterwine_data;
    case 1  %'rs'  
        dataOut_dc = LSY_RSDec(deinterwine_data,nn,kk);
%         [num,Ber(i)]=biterr(dataOut_dc,dataIn);
    case 2  %'juanji'
        dataOut_dc = vitdec(deinterwine_data,trel,tbdepth,'trunc','hard');             %Viterbi译码 
        %cont表示连续调用模式(此时在第一个被解码的符号出现在输出中之前，会出现一个等于输入tbdepth符号的延迟。误比特计算要用下面这种方式) trunc表示截断模式
end
%% 检错解码
switch CRC_CODE
    case 0
        dataOut = dataOut_dc;
    case 1  %(237,225)检错编码
        [dataOut_crc, dataOut_crc_err] = crcdetector(dataOut_dc');%输入为列向量 dataOut_crc CRC解码输出 dataOut_crc_err错误个数
        CRC(i) = sum(dataOut_crc_err);
        dataOut = dataOut_crc';
end
% [num,Ber(i)] = biterr(dataOut,dataIn); %计算误比特率

        % 计算误bit数(包括同步序列和消息序列，不包括冗余序列)
        errBitNum = sum(dataIn ~= dataOut);
        sumErrBit = sumErrBit + errBitNum;
        
        % 计算误包率
        orignalPacket = dataIn(SYNC_BIT_NUM+1 : end);
        dmdPacket = dataOut(SYNC_BIT_NUM+1 : end);        
        orignalPacketMat = reshape(orignalPacket , MSG_BIT_NUM/PACKET_NUM , PACKET_NUM);    % bit序列转换成包，每一列一包
        dmdPacketMat = reshape(dmdPacket , MSG_BIT_NUM/PACKET_NUM , PACKET_NUM);            % bit序列转换成包，每一列一包
        % 统计误包(比较一包里面的bit是否与发送端向相同，不同加1，相同加0)
        sumErrPacket = ...
            (sum(orignalPacketMat(:, 1) ~= dmdPacketMat(:, 1))>0)+...
            (sum(orignalPacketMat(:, 2) ~= dmdPacketMat(:, 2))>0)+...
            (sum(orignalPacketMat(:, 3) ~= dmdPacketMat(:, 3))>0)+sumErrPacket;
        
        [ii jj]
    end 
    % 计算误bit率
    ber(ii) = sumErrBit / (frameNum*TX_BIT_NUM);
    berPacket(ii) = sumErrPacket/(frameNum*PACKET_NUM);
end

%% 误码输出
figure
semilogy(snr,ber,'-*')
save MSK+RS+CCSK+FH ber
hold on
grid on

BER_qpsk = duibizu_qpsk(dataIn,snr,PACKET_NUM);
save MSK+RS+CCSK BER_qpsk
semilogy(snr,BER_qpsk,'-o')
hold on

BER_16qam = duibizu_16qam(dataIn,snr,PACKET_NUM);
save MSK+CCSK BER_16qam
semilogy(snr,BER_16qam,'-x')
hold on

BER_64qam = duibizu_64qam(dataIn,snr,PACKET_NUM);
save MSK BER_64qam
semilogy(snr,BER_64qam,'-s')
hold on

grid on
xlabel('Eb/No (dB)')%如果不进行转换就是SNR
ylabel('误比特率')%如果不进行转换就是误码率

legend('MSK+RS+CCSK+FH','MSK+RS+CCSK','MSK+CCSK','MSK')

toc