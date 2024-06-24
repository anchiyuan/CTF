clc; clear;
close all;

%% RIR parameter %%
SorNum = 2;                                              % source number
MicNum = 30;                                             % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)
Ts = 1/fs;                                               % Sample period (s)

% ULA %
MicStart = [1, 1.5, 1];
spacing = 0.02;
MicPos = zeros(MicNum, 3);
for i = 1:MicNum
    MicPos(i, :) = [MicStart(1, 1)+(i-1)*spacing MicStart(1, 2) MicStart(1, 3)];
end

SorPos = [3.5, 2.6, 1 ; 0.5, 4, 1];                        % source position (m)
room_dim = [5, 6, 2.5];                                  % Room dimensions [x y z] (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reverberation_time = 0.2;                                % Reverberation time (s)
points_rir = 4096;                                       % Number of rir points (需比 reverberation time 還長)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 1;                                           % Disable high-pass filter

% 畫空間圖 %
figure(1);
plot3( [0 room_dim(1, 1) room_dim(1, 1) 0 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1)], ...
       [0 0 room_dim(1, 2) room_dim(1, 2) 0 0 0 room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) 0 0 0], ...
       [0 0 0 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0] , 'k')
hold on
plot3(MicPos(:, 1), MicPos(:, 2), MicPos(:, 3), 'r.', 'MarkerSize', 10)
hold on
plot3(SorPos(:, 1), SorPos(:, 2), SorPos(:, 3), '*', 'MarkerSize', 20)
hold off
xlabel('x\_axis')
ylabel('y\_axis')
zlabel('z\_axis')
title('空間圖')
shg

%% generate ground-truth RIR (h) %%
% 產生 RIR 和存.mat 檔 %
% h = zeros(MicNum, SorNum, points_rir);
% h(:, 1, :) = rir_generator(c, fs, MicPos, SorPos(1, :), room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
% h(:, 2, :) = rir_generator(c, fs, MicPos, SorPos(2, :), room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
% rir_filename_str = ['h_TIKR\h_', string(reverberation_time), 'x', string(MicNum), 'x', string(SorNum), 'x', string(points_rir), '.mat'];
% rir_filemane = join(rir_filename_str, '');
% save(rir_filemane, 'h')

% load RIR 的 .mat 檔 %
rir_filename_str = ['h_TIKR\h_', string(reverberation_time), 'x', string(MicNum), 'x', string(SorNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

look_mic = 10;
% 畫 ground-truth RIR time plot %
figure(2)
subplot(2, 1, 1);
plot(squeeze(h(look_mic, 1, :)).', 'r');
title('ground-truth RIR')
xlabel('points')
ylabel('amplitude')

subplot(2, 1, 2);
plot(squeeze(h(look_mic, 2, :)).', 'r');
title('ground-truth RIR')
xlabel('points')
ylabel('amplitude')
shg

%% window parameter %%
NFFT = points_rir;
hopsize = NFFT/4;
win = hamming(NFFT);
frequency = NFFT/2 + 1;
freqs_vector = linspace(0, fs/2, frequency);

%% RIR 轉頻域 (H) %%
H = zeros(MicNum, SorNum, frequency);
H_temp = fft(squeeze(h(:, 1, :)), NFFT, 2);
H(:, 1, :) = H_temp(:, 1:frequency);
H_temp = fft(squeeze(h(:, 2, :)), NFFT, 2);
H(:, 2, :) = H_temp(:, 1:frequency);

%% 讀音檔 or 產生 white noise source (source) %%
Second = 23;
SorLen =  Second*fs;

% load speech source %
source_transpose = zeros(SorLen, SorNum);
[source_transpose(:, 1), ~] = audioread('245.wav', [1, SorLen]);
[source_transpose(:, 2), ~] = audioread('313.wav', [1, SorLen]);% speech source
source = source_transpose.';

% load white noise source %
% source = wgn(1, SorLen, 0);                                    % white noise source

%% compute source signal for frequency (S) %%
% source 轉頻域 %
source_transpose = source.';
[S, ~, ~] = stft(source_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

NumOfFrame = size(S, 2);
NumOfFrame_vector = 1:1:NumOfFrame;

%% 產生麥克風訊號，先在時域上 convolution 再做 stft (y Y) %%
% convolution source and RIR %
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(squeeze(h(i, 1, :)).', source(1, :));
    as(i, :) = as(i, :) + conv(squeeze(h(i, 2, :)).', source(2, :));
end

y =as;

% y 轉頻域 %
y_transpose = y.';
[Y, ~, ~] = stft(y_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% TIKR %%
beta = 1e-7;
S_TIKR = zeros(frequency, NumOfFrame, SorNum);
for n = 1:frequency
    for FrameNo = 1:NumOfFrame
        S_TIKR(n, FrameNo, :) = inv(H(:, :, n)'*H(:, :, n) + beta*eye(SorNum))*H(:, :, n)'*squeeze(Y(n, FrameNo, :));
    end
end

source_TIKR_transpose = istft(S_TIKR, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, ConjugateSymmetric=true, FrequencyRange='onesided');
source_TIKR = source_TIKR_transpose.';

wav_filename_str = ['wav_TIKR\source2_', string(beta), '.wav'];
wav_filename = join(wav_filename_str, '');
audiowrite(wav_filename, source_TIKR(2, :), fs)
