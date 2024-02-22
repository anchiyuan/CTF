clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

%% RIR parameter %%
SorNum = 1;                                              % source number
MicNum = 6;                                             % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)
Ts = 1/fs;                                               % Sample period (s)
points_rir = 2048;                                      % Number of rir points (需比 reverberation time 還長)

%% compute ground-truth CTF (H) %%
NFFT = 1024;
hopsize = 256;

% windows %
win = hamming(NFFT);
osfac = round(NFFT/hopsize);

% compute H %
frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);

L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);


%% 讀音檔 or 產生 white noise source (source) %%
Second = 10;
SorLen =  Second*fs;

% load source %
[source_transpose, fs] = audioread('exp_wav\loudspeasker_1.wav', [1, SorLen]);    % speech source
source = source_transpose.';

% source = wgn(1, SorLen, 0);                                    % white noise source

%% compute source signal for frequency (S) %%
% source 轉頻域 %
source_transpose = source.';
[S, ~, S_t_vector] = stft(source_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

NumOfFrame = size(S_t_vector, 1);
NumOfFrame_vector = 1:1:NumOfFrame;
S_dB = mag2db(abs(S));
% 畫 source frequency plot %
% figure(4);
% mesh(NumOfFrame_vector, freqs_vector, S_dB)
% colorbar
% view(2)
% xlim([1 NumOfFrame])
% title('S')
% xlabel('frame')
% ylabel('frequency(Hz)')
% shg

%% RIR mix source 先在時域上 convolution 再做 stft (y and Y) %%
for i = 1:MicNum
    y_nodelay_str = ['exp_wav\', string(i), 'th.wav'];
    y_nodelay_filename = join(y_nodelay_str, '');
    [y_nodelay(i, :), fs] = audioread( y_nodelay_filename, [1, SorLen]); 
end
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, SorLen);
y_delay(:, extra_delay_y+1:end) = y_nodelay(:, 1:SorLen-extra_delay_y);


% y 轉頻域 %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');


%% predict CTF with Kalman filter (A) %%
start_ini_frame = L;    % wpe 第12秒開始穩定
ini_frame = NumOfFrame;

% Kalman filter %
A = zeros(MicNum, L, frequency);
tic
parfor i = 1:MicNum
    for n =1:frequency
        weight = zeros(L, 1);
        P = 0.5*eye(L);    % error covariance matrix %
        K = zeros(L, 1);    % Kalman gain %
        R = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_ini_frame:ini_frame
            % no time update only have measurement update %
            K = P*flip(S(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(S(n, FrameNo-L+1:FrameNo)))*P*flip(S(n, FrameNo-L+1:FrameNo).') + R);
            weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(S(n, FrameNo-L+1:FrameNo)))*weight);
            P = P - K*conj(flip(S(n, FrameNo-L+1:FrameNo)))*P;
        end
    
        A(i, :, n) = weight';
    end

end
toc

look_mic = 3;
A_dB = mag2db(abs(A));
% 畫 A frequency plot %
figure(1);
mesh(L_vector, freqs_vector, squeeze(A_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% TF estimate %%
TF =  tfestimate(source, y_nodelay.', hamming(points_rir), points_rir*0.75, points_rir);
TF = TF.';
TF = [TF flip(conj(TF(:, 2:end-1)))];
tf = ifft(TF, points_rir, 2,'symmetric');

figure(2)
plot(tf(look_mic, :));
xlabel('points')
ylabel('amplitude')
shg

%% A 轉回時域 (A_tdomain) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
ratio_A_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_tdomain(i, :) = max(abs(tf(i, :)))/max(abs(A_tdomain(i, hopsize*(osfac-1)+1:end)));
end

A_tdomain = A_tdomain.*ratio_A_tdomain;

A_tdomain_forplot = A_tdomain(look_mic, hopsize*(osfac-1)+1:end);
% 畫 A_tdomain time plot
figure(7)
plot(tf(look_mic, :), 'r');
hold on
plot(A_tdomain_forplot, 'b');
hold off
legend('MTF', 'CTF')
xlabel('points')
ylabel('amplitude')
shg

%%  dereverb with MINT (source_MINT) %%
A_tdomain_forMINT = A_tdomain(:, hopsize*(osfac-1)+1:end);
g_len = 6000;
weight_len = floor(g_len/4);
dia_load_MINT = 10^(-7);
source_MINT = MINT(A_tdomain_forMINT, y_nodelay, g_len, weight_len, dia_load_MINT);

% adjust source_predict 的最高點使之與 source 的一樣 %
source_max  = max(abs(source(1, :)));
source_MINT_max  = max(abs(source_MINT(1, :)));
ratio_source_MINT = source_max/source_MINT_max;
source_MINT = source_MINT.*ratio_source_MINT;

% 畫 source_predict time plot %
figure(8)
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
% save all wav %
% audiowrite('wav\source_all.wav', source(1, :), fs)
% 
% ratio_y_nodelay = 0.8 / max(abs(y_nodelay(look_mic, :))) ;
% y_filemane_str = ['wav\y_nodelay_all-', string(reverberation_time), '.wav'];
% y_filemane = join(y_filemane_str, '');
% audiowrite(y_filemane, y_nodelay(look_mic, :)*ratio_y_nodelay, fs)
% 
% source_MINT_filemane_str = ['wav\source_predict_all_Kalman_MINT-', string(reverberation_time), '.wav'];
% source_MINT_filemane = join(source_MINT_filemane_str, '');
% audiowrite(source_MINT_filemane, source_MINT(1, :), fs)

% save partial wav %
point_start_save = 18*fs;

audiowrite('wav\source_partial.wav', source(1, point_start_save:end), fs)

ratio_y_nodelay = 0.8 / max(abs(y_nodelay(look_mic, point_start_save:end))) ;
y_filemane_str = ['wav\y_nodelay_partial-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_nodelay(look_mic, point_start_save:end)*ratio_y_nodelay, fs)

ratio_y_wpe = 0.8 / max(abs(y_wpe(look_mic, point_start_save:end))) ;
y_filemane_str = ['wav\y_wpe_partial-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_wpe(look_mic, point_start_save:end)*ratio_y_wpe, fs)

source_MINT_filemane_str = ['wav\source_predict_partial_Kalman_MINT-', string(reverberation_time), '.wav'];
source_MINT_filemane = join(source_MINT_filemane_str, '');
audiowrite(source_MINT_filemane, source_MINT(1, point_start_save:end), fs)

fprintf('done\n')
