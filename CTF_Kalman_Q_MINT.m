clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

%% RIR parameter %%
SorNum = 1;                                              % source number
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

SorPos = [2, 2.6, 1];                                    % source position (m)
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

%% generate ground-truth RIR (h) %%
% 產生 RIR 和存.mat 檔 %
% h = rir_generator(c, fs, MicPos, SorPos, room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
% rir_filename_str = ['h\h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
% rir_filemane = join(rir_filename_str, '');
% save(rir_filemane, 'h')

% load RIR 的 .mat 檔 %
rir_filename_str = ['h\h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

look_mic = 10;
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
% 畫 ground-truth RIR time plot %
figure(1)
plot(h(look_mic, :), 'r');
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('ground-truth RIR')
xlabel('points')
ylabel('amplitude')
shg

%% window parameter %%
NFFT = 1024;
hopsize = 256;

% windows %
win = hamming(NFFT);
osfac = round(NFFT/hopsize);

frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);    % (len(win) + len(win) - 1) + points_rir - 1
L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);

%% 讀音檔 or 產生 white noise source (source) %%
Second = 23;
SorLen =  Second*fs;

% load speech source %
[source_transpose, fs] = audioread('245.wav', [1, SorLen]);    % speech source
source = source_transpose.';

% load white noise source %
% source = wgn(1, SorLen, 0);                                    % white noise source

%% compute source signal for frequency (S) %%
% source 轉頻域 %
source_transpose = source.';
[S, ~, S_t_vector] = stft(source_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

NumOfFrame = size(S_t_vector, 1);
NumOfFrame_vector = 1:1:NumOfFrame;

%% 產生麥克風訊號，先在時域上 convolution 再做 stft (y_nodelay y_delay Y_delay ) %%
% convolution source and RIR %
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(h(i, :), source);
end

extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, SorLen);
y_delay(:, extra_delay_y+1:end) = as(:, 1:SorLen-extra_delay_y);
y_nodelay = as(:, 1:SorLen);

% y_delay 轉頻域 %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPE (y_wpe) %%
% % do wpe %
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
% y_wpe = y_wpe.';
% 
% % 存 wpe mat %
% y_wpe_filename_str = ['y\y_wpe_', string(reverberation_time), '.mat'];
% y_wpe_filename = join(y_wpe_filename_str, '');
% save(y_wpe_filename, 'y_wpe')

% load y_wpe %
y_wpe_filename_str = ['y\y_wpe_', string(reverberation_time), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

%% DAS beamformer (Y_DAS) %%
% 算 mic 與 source 之距離 %
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((SorPos - MicPos(i, :)).^2));
end

% 算 freefield ATF %
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

%% predict CTF with Kalman filter (exist process noise) (A) %%
start_ini_frame = L;
ini_frame = NumOfFrame;

% Kalman filter %
A = zeros(MicNum, L, frequency);
tic
parfor i = 1:MicNum
    for n =1:frequency
        weight = zeros(L, 1);
        P = 0.5*eye(L);    % error covariance matrix %
        K = zeros(L, 1);    % Kalman gain %
        Q = 10^(-4)*eye(L);    % process noise covariance matrix %
        R = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_ini_frame:ini_frame
            % time update %
            P = P + Q;

            % measurement update %
            K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
            weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
            P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;

        end

        A(i, :, n) = weight';
    end

end
toc

%% A 轉回時域且算 NRMSPM (A_tdomain NRMSPM) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
A_tdomain = A_tdomain(:, hopsize*(osfac-1)+1:end);
ratio_A_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_tdomain(i, :) = max(abs(h(i, :)))/max(abs(A_tdomain(i, :)));
end

A_tdomain = A_tdomain.*ratio_A_tdomain;

% 畫 A_tdomain time plot
figure(2)
plot(h(look_mic, :), 'r');
hold on
plot(A_tdomain(look_mic, :), 'b');
hold off
title('estimated RIR')
legend('ground-truth RIR', 'estimated RIR')
xlabel('points')
ylabel('amplitude')
shg

%  算 NRMSPM %
h_NRMSPM = reshape(h.', [MicNum*points_rir 1]);
aa_NRMSPM = reshape(A_tdomain.', [MicNum*points_rir 1]);
NRMSPM = 20*log10(norm(h_NRMSPM-h_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(h_NRMSPM));

%%  dereverb with MINT (source_MINT) %%
g_len = floor(points_rir/1000)/4*3000;
weight_len = floor(g_len/4);
dia_load_MINT = 10^(-7);
source_MINT = MINT(A_tdomain, y_nodelay, g_len, weight_len, dia_load_MINT);

% adjust source_predict 的最高點使之與 source 的一樣 %
source_max  = max(abs(source(1, :)));
source_MINT_max  = max(abs(source_MINT(1, :)));
ratio_source_MINT = source_max/source_MINT_max;
source_MINT = source_MINT.*ratio_source_MINT;

% 畫 source_predict time plot %
figure(3)
plot(source(1, :), 'r');
hold on
plot(source_MINT(1, :), 'b');
hold off
title('source\_MINT')
xlabel('points')
ylabel('magnitude')
legend('source', 'source\_MINT')
shg

%% save .wav 檔 %%
% save partial wav %
point_start_save = 18*fs;

audiowrite('wav\source_partial.wav', source(1, point_start_save:end), fs)

ratio_y_nodelay = 0.8 / max(abs(y_nodelay(look_mic, point_start_save:end))) ;
y_filemane_str = ['wav\y_nodelay_partial_', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_nodelay(look_mic, point_start_save:end)*ratio_y_nodelay, fs)

ratio_y_wpe = 0.8 / max(abs(y_wpe(look_mic, point_start_save:end))) ;
y_filemane_str = ['wav\y_wpe_partial_', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_wpe(look_mic, point_start_save:end)*ratio_y_wpe, fs)

source_MINT_filemane_str = ['wav\source_predict_partial_Kalman_Q_MINT_', string(reverberation_time), '.wav'];
source_MINT_filemane = join(source_MINT_filemane_str, '');
audiowrite(source_MINT_filemane, source_MINT(1, point_start_save:end), fs)

fprintf('done\n')
