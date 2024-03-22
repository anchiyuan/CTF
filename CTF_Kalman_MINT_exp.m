clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

%% RIR parameter %%
SorNum = 1;
c = 343;
fs = 48000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MicNum = 6;           % 實驗麥克風數量
SpeakerNum = 6;       % 實驗喇叭數量
look_speaker = 3;     % source 為哪一個喇叭
look_mic = 1;         % 指定想要畫圖之麥克風
points_rir = 8192;    % 自行設定想要輸出的 RIR 長度
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
Second_speaker = 15;    % 每顆喇叭播放時間 (10+5)
Second = 8;             % 使用時間長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*fs;

% load source %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[source_transpose, Fs] = audioread('wav_exp\loudspeaker_1.wav', [1, SorLen]);    % speech source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source = source_transpose.';

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
    y_nodelay_str = ['wav_exp\', string(i-1), 'th.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_nodelay_filename = join(y_nodelay_str, '');
    [y_nodelay(i, :), Fs] = audioread( y_nodelay_filename, [(look_speaker-1)*Second_speaker*fs+1, (look_speaker-1)*Second_speaker*fs+SorLen]);
end

% delay y_nodelay to get y_delay %
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, SorLen);
y_delay(:, extra_delay_y+1:end) = y_nodelay(:, 1:SorLen-extra_delay_y);

% y_delay 轉頻域 to get Y_delay %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPE (y_wpe) %%
% % do wpe %
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
% y_wpe = y_wpe.';
% 
% % 存 wpe mat %
% y_wpe_str = ['y_exp\y_wpe_', string(look_speaker), '.mat'];
% y_wpe_filename = join(y_wpe_str, '');
% save(y_wpe_filename, 'y_wpe')

% load y_wpe %
y_wpe_str = ['y_exp\y_wpe_', string(look_speaker), '.mat'];
y_wpe_filename = join(y_wpe_str, '');
load(y_wpe_filename);

%% DAS beamformer (Y_DAS) %%
% mic 與 source 之距離 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if look_speaker == 1
    distance = [1.1; 1.153; 1.203; 1.253; 1.303; 1.36];    % 第一顆喇叭
else
    distance = [0.828; 0.83; 0.834; 0.84; 0.852; 0.87];    % 第三顆喇叭
end
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

%% A 轉回時域且算 NRMSPM (A_tdomain NRMSPM) %%
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

% 畫 A_tdomain time plot
figure(2)
plot(tf(look_mic, :), 'r');
hold on
plot(-A_tdomain(look_mic, :), 'b');
hold off
legend('MTF', 'CTF')
xlabel('points')
ylabel('amplitude')
shg

% 算 NRMSPM %
h_NRMSPM = reshape(tf.', [MicNum*points_rir 1]);
aa_NRMSPM = reshape(A_tdomain.', [MicNum*points_rir 1]);
NRMSPM = 20*log10(norm(h_NRMSPM-h_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(h_NRMSPM));

%% 檢查 direct sound 有無 match %%
[~, argmax_h] = max(abs(tf.'));
[~, argmax_A_tdomain] = max(abs(A_tdomain.'));

NRMSPM_in = zeros(MicNum, 1);
for i = 1:MicNum
    NRMSPM_in(i, :) = 20*log10(norm(tf(i, :).'-tf(i, :)*A_tdomain(i, :).'/(A_tdomain(i, :)*A_tdomain(i, :).')*A_tdomain(i, :).')/norm(tf(i, :).'));
end

%% 檢查 ATF 有無 match %%
ATF = fft(tf, points_rir, 2);
ATF_estimated = fft(A_tdomain, points_rir, 2);

figure(3)
semilogx(linspace(0, fs/2, points_rir/2+1), abs(ATF(look_mic, 1:points_rir/2+1)), 'r');
hold on
semilogx(linspace(0, fs/2, points_rir/2+1), abs(ATF_estimated(look_mic, 1:points_rir/2+1)), 'b');
hold off
legend('MTF', 'CTF')
xlabel('frequency')
ylabel('magnitude')
shg

fprintf('done\n')
