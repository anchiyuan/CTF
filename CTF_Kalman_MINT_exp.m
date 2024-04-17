clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

%% RIR parameter %%
SorNum = 1;
c = 343;
Fs = 48000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 16000;           % 欲 resample 成的取樣頻率
MicNum = 6;           % 實驗麥克風數量
look_mic = 1;         % 指定想要畫圖之麥克風
points_rir = 2048;    % 自行設定想要輸出的 RIR 長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% window parameter %%
NFFT = 1024;
hopsize = 256;

win = hamming(NFFT);
osfac = round(NFFT/hopsize);

frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);
L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);

%% read source 音檔 (source) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Second = 28;             % 使用時間長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*Fs;
sorLen =  Second*fs;

% load source %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[source_transpose, ~] = audioread('wav_exp\white_noise_30s.wav', [1, SorLen]);    % speech source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source = source_transpose.';

% resample %
source = resample(source, 1, Fs/fs);

%% source 做 stft (S) %%
source_transpose = source.';
[S, ~, ~] = stft(source_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

NumOfFrame = size(S, 2);
NumOfFrame_vector = 1:1:NumOfFrame;

%% read mic 音檔再做 stft (y_nodelay, y_delay and Y_delay) %%
% load y_nodelay %
y_nodelay = zeros(MicNum, SorLen);
for i = 1:MicNum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_nodelay_str = ['wav_exp\', string(i), '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_nodelay_filename = join(y_nodelay_str, '');
    [y_nodelay(i, :), ~] = audioread( y_nodelay_filename, [1, SorLen]);
end

% resample %
y_nodelay = resample(y_nodelay, 1, Fs/fs, Dimension=2);

% delay y_nodelay to get y_delay %
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, sorLen);
y_delay(:, extra_delay_y+1:end) = y_nodelay(:, 1:sorLen-extra_delay_y);

% y_delay 轉頻域 to get Y_delay %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPE (y_wpe) %%
% % do wpe %
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
% y_wpe = y_wpe.';
% 
% % 存 wpe mat %
% y_wpe_str = ['y_exp\y_wpe_', string(fs),'.mat'];
% y_wpe_filename = join(y_wpe_str, '');
% save(y_wpe_filename, 'y_wpe')

% load y_wpe %
y_wpe_str = ['y_exp\y_wpe_', string(fs),'.mat'];
y_wpe_filename = join(y_wpe_str, '');
load(y_wpe_filename);

%% DAS beamformer (Y_DAS) %%
% mic 與 source 之距離 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distance = [0.83; 0.832; 0.832; 0.84; 0.852; 0.874];    % 第三顆喇叭

% SorPos = [1.4126*sind(-6.8484), 1.4126*cosd(-6.8484), 0];
% distance = zeros(MicNum, 1);
% for i = 1:MicNum
%     distance(i, :) = norm(SorPos - [(i-1)*0.07, 0, 0]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 算 a %
a = zeros(MicNum, SorNum, frequency);
for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, :, n) = exp(-1j*omega/c*distance)./distance;
end

% 算 DAS weight %
w = a/MicNum;

% 算 Y_DAS %
y_wpe_transpose = y_wpe.';
[Y_wpe, ~, ~] = stft(y_wpe_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

Y_DAS = zeros(frequency, NumOfFrame);
for FrameNo= 1:NumOfFrame
    for n = 1: frequency
         Y_DAS(n, FrameNo) = w(:, :, n)'*squeeze(Y_wpe(n, FrameNo, :));
    end  

end

%% predict CTF with Kalman filter (A) %%
start_ini_frame = L;
ini_frame = NumOfFrame;
A = zeros(MicNum, L, frequency);

% Kalman stationary filter %
tic
parfor i = 1:MicNum
    for n =1:frequency
        weight = zeros(L, 1);
        P = 0.5*eye(L);    % error covariance matrix %
        K = zeros(L, 1);    % Kalman gain %
        R = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_ini_frame:ini_frame
            % no time update only have measurement update %
            K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
            weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
            P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;
        end
    
        A(i, :, n) = weight';
    end

end
toc

% Kalman nonstationary filter %
% tic
% parfor i = 1:MicNum
%     for n =1:frequency
%         weight = zeros(L, 1);
%         P = 0.5*eye(L);    % error covariance matrix %
%         K = zeros(L, 1);    % Kalman gain %
%         Q = 10^(-4)*eye(L);    % process noise covariance matrix %
%         R = 10^(-3);    % measurement noise covariance matrix %
%         for FrameNo = start_ini_frame:ini_frame
%             % time update %
%             P = P + Q;
% 
%             % measurement update %
%             K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
%             weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
%             P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;
% 
%         end
% 
%         A(i, :, n) = weight';
%     end
% 
% end
% toc

%% TF estimate (tf and t60) %%
TF =  tfestimate(source, y_nodelay.', hamming(points_rir), points_rir*0.75, points_rir);
TF = TF.';
TF = [TF flip(conj(TF(:, 2:end-1)))];
tf = ifft(TF, points_rir, 2,'symmetric');

e_total = sum(tf.^2, 2);
t60 = zeros(MicNum, 1);
for i = 1:MicNum
    edc = zeros(points_rir, 1);
    for j = 1:points_rir
        edc(j, :) = 10*log10(sum(tf(i, j:end).^2)/e_total(i, :));
    end
    
    edc = edc + 60;
    [~, t60_point] = min(abs(edc));
    t60(i, :) = t60_point/fs;
end

figure(1)
plot(tf(look_mic, :));
title('tfestimate')
xlabel('points')
ylabel('amplitude')
shg

%% A 轉回時域畫圖 (A_tdomain) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
A_tdomain = A_tdomain(:, hopsize*(osfac-1)+1:end);
ratio_A_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_tdomain(i, :) = max(abs(tf(i, :)))/max(abs(A_tdomain(i, :)));
end

A_tdomain = A_tdomain.*ratio_A_tdomain;

look_mic = 1;
% 畫 A_tdomain time plot
figure(2)
plot(tf(look_mic, :), 'r');
hold on
plot(-A_tdomain(look_mic, :), 'b');
hold off
legend('tfestimate', 'CTF')
title('RIR')
xlabel('points')
ylabel('amplitude')
shg

%% 算 NRMSPM %%
tf_NRMSPM = reshape(tf.', [MicNum*points_rir 1]);
A_tdomain_NRMSPM = reshape(A_tdomain.', [MicNum*points_rir 1]);
NRMSPM = 20*log10(norm(tf_NRMSPM-tf_NRMSPM.'*A_tdomain_NRMSPM/(A_tdomain_NRMSPM.'*A_tdomain_NRMSPM)*A_tdomain_NRMSPM)/norm(tf_NRMSPM));

NRMSPM_in = zeros(MicNum, 1);
for i = 1:MicNum
    NRMSPM_in(i, :) = 20*log10(norm(tf(i, :).'-tf(i, :)*A_tdomain(i, :).'/(A_tdomain(i, :)*A_tdomain(i, :).')*A_tdomain(i, :).')/norm(tf(i, :).'));
end

%% 檢查 direct sound 有無 match %%
[~, argmax_tf] = max(abs(tf.'));
[~, argmax_A_tdomain] = max(abs(A_tdomain.'));

%% 檢查 ATF 有無 match %%
ATF = fft(tf, points_rir, 2);
ATF_estimated = fft(A_tdomain, points_rir, 2);

figure(3)
semilogx(linspace(0, fs/2, points_rir/2+1), abs(ATF(look_mic, 1:points_rir/2+1)), 'r');
hold on
semilogx(linspace(0, fs/2, points_rir/2+1), abs(ATF_estimated(look_mic, 1:points_rir/2+1)), 'b');
hold off
legend('tfestimate', 'CTF')
title('ATF')
xlabel('frequency')
ylabel('magnitude')
shg

fprintf('done\n')
