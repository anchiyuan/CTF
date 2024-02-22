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
reverberation_time = 0.4;                                % Reverberation time (s)
points_rir = 8192;                                       % Number of rir points (需比 reverberation time 還長)
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 0;                                           % Disable high-pass filter

%% generate ground-truth RIR (h) %%
% load RIR 的 .mat 檔 %
rir_filename_str = ['h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

look_mic = 10;
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
% 畫 ground-truth RIR time plot %
figure(1)
plot(h(look_mic, :), 'r');
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h')
xlabel('points')
ylabel('amplitude')
shg

%% compute ground-truth CTF (H) %%
NFFT = 1024;
hopsize = 256;

% windows %
win = hamming(NFFT);
osfac = round(NFFT/hopsize);

% compute H %
frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);
H = zeros(frequency, L, MicNum);
for k = 1:frequency
    zeta = xcorr(win,win).*exp(1i*2*pi*(k-1)*(-NFFT+1:NFFT-1)'/NFFT);
    for i = 1:MicNum
        H_temp = conv(h(i, :), zeta.');
        H(k, :, i) = H_temp(:, hopsize:hopsize:end);
    end
end

%% reconstruct RIR from H (h_reconstruct) %%
h_reconstruct = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, H);

% 算 reconstruct ratio %
ratio_h_reconstruct = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_h_reconstruct(i, :) = max(abs(h(i, :)))/max(abs(h_reconstruct(i, :)));
end

% h_reconstruct_forplot = h_reconstruct(:, hopsize*(osfac-1)+1:end);
% h_reconstruct_forplot = h_reconstruct_forplot.*ratio_h_reconstruct;
% 畫 reconstructed RIR time plot %
% figure(2)
% plot(h(look_mic, :), 'r');
% hold on
% plot(h_reconstruct_forplot(look_mic, :), 'b');
% hold off
% ylim([h_yaxis_underlimit h_yaxis_upperlimit])
% title('reconstruct h from H')
% legend('h', 'h\_reconstruct')
% xlabel('points')
% ylabel('amplitude')
% shg

% 把 H 調成由他 reconstruct 之 h 會有正確大小後再算他的 PSD %
PSD_H = zeros(MicNum, 1);
for i = 1:MicNum
    H(:, :, i) = H(:, :, i)*ratio_h_reconstruct(i, :);
    PSD_H(i, :) = sum(abs(H(:, :, i)).^2, "all");
end

L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);
H_dB = mag2db(abs(H));
% 畫 ground-truth RIR frequency plot %
figure(3);
mesh(L_vector, freqs_vector, H_dB(:, :, look_mic))
colorbar
view(2)
xlim([1 L])
title('H')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% 讀音檔 or 產生 white noise source (source) %%
Second = 23;
SorLen =  Second*fs;

% load source %
[source_transpose, fs] = audioread('245.wav', [1, SorLen]);    % speech source
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
% convolution source and RIR %
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(h(i, :), source);
end

extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, SorLen);
y_delay(:, extra_delay_y+1:end) = as(:, 1:SorLen-extra_delay_y);
y_nodelay = as(:, 1:SorLen);

% y 轉頻域 %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

y_nodelay_transpose = y_nodelay.';
[Y_nodelay, ~, ~] = stft(y_nodelay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPD and RLS %%
b = 2;
frequency_range = [800, 1500, 3000];
Lw = [25, 20, 15, 10];    % 0–800, 800–1500, 1500–3000, and 3000–4000 Hz.
Lf = 2;    %    length of WPD for compute time-varying power of the source signal
diaload_R_bar = 1e-4;
lambda = 0.99;
start_ini_frame = Lw(1) + L;
ini_frame = NumOfFrame;
reference_mic = 1;
S_WPD = Y_nodelay(:, :, reference_mic);
S_WPD = [S_WPD zeros(frequency, floor(Lf/2))];    %  sigma 如果有考慮到 future 時才不會 error
A = zeros(MicNum, L, frequency);
error = zeros(MicNum, (ini_frame-start_ini_frame+1), frequency);

% compute distance for a %
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((SorPos - MicPos(i, :)).^2));
end

% iteration process %
for n = 1: frequency
    % choose Lw for different frequency %
    if (freqs_vector(n) <= frequency_range(1))
        Lw_choose = Lw(1);
    elseif (freqs_vector(n) <= frequency_range(2))
        Lw_choose = Lw(2);
    elseif (freqs_vector(n) <= frequency_range(3))
        Lw_choose = Lw(3);
    else
        Lw_choose = Lw(4);
    end
    
    % initialize R_bar %
    R_bar = zeros(MicNum*(Lw_choose-b+2), MicNum*(Lw_choose-b+2));

    % construct a_bar %
    omega = 2*pi*freqs_vector(n);
    a = exp(-1j*omega/c*distance)./distance;
    a_bar = zeros(MicNum*(Lw_choose-b+2), 1);
    a_bar(1:MicNum, :) = a;
    
    % initialize RLS variable %
    weight = zeros(L, MicNum);
    P = (1e-2)^-1*eye(L);

    for FrameNo = start_ini_frame:ini_frame
        % cascade right now frame mic signal %
        Y_bar_now = squeeze(Y_nodelay(n, FrameNo, :));
        for frame = b:Lw_choose
            Y_bar_now = cat(1, Y_bar_now, squeeze(Y_nodelay(n, FrameNo-frame, :)));
        end

        % compute sigma %
%         sigma = mean(abs(S_WPD(n, FrameNo-floor(Lf/2):FrameNo+floor(Lf/2))).^2, 'all');    % use past and now and future to compute sigma
        sigma = mean(abs(S_WPD(n, FrameNo-Lf+1:FrameNo)).^2, 'all');    % use only past and now to compute sigma

        % compute WPD weight %
        R_bar = R_bar + Y_bar_now*Y_bar_now'/sigma;
        w_bar = inv(R_bar + diaload_R_bar*eye(MicNum*(Lw_choose-b+2)))*a_bar/(a_bar'*inv(R_bar + diaload_R_bar*eye(MicNum*(Lw_choose-b+2)))*a_bar);

        % use WPD compute source siganal %
        for l = 1:L
            Y_bar = squeeze(Y_nodelay(n, FrameNo-l+1, :));
            for frame = b:Lw_choose
                Y_bar = cat(1, Y_bar, squeeze(Y_nodelay(n, FrameNo-l+1-frame, :)));
            end

            S_WPD(n, FrameNo-l+1) = w_bar'*Y_bar;
        end

        % RLS %
        S_WPD_choose = flip(S_WPD(n, FrameNo-L+1:FrameNo).');
        Y_delay_choose = squeeze(Y_delay(n, FrameNo, :));
        error(:, FrameNo-start_ini_frame+1, n) = Y_delay_choose' - S_WPD_choose'*weight;
        k_RLS = lambda^-1*P*S_WPD_choose/(1 + lambda^-1*S_WPD_choose'*P*S_WPD_choose);
        weight = weight + k_RLS*error(:, FrameNo-start_ini_frame+1, n).';
        P = lambda^-1*P - lambda^-1*k_RLS*S_WPD_choose'*P;

    end

    A(:, :, n) = weight';
end

error_abs_mean = mean(mean(abs(error), 1), 3);
% 畫 error curve %
figure(5)
plot(error_abs_mean)
title('mean absolute error')
xlabel('iteration')
ylabel('magnitude')

% 調整 A 的 PSD 使之與 H 的一樣 %
for i = 1:MicNum
    ratio_PSD = PSD_H(i, :)/sum(abs(A(i, :, :)).^2, "all");
    A(i, : ,:) = A(i, : ,:)*sqrt(ratio_PSD);
end

A_dB = mag2db(abs(A));
% 畫 A frequency plot %
figure(6);
mesh(L_vector, freqs_vector, squeeze(A_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A\_ini')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% A 轉回時域 (A_tdomain) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);

A_tdomain_forplot = A_tdomain(look_mic, hopsize*(osfac-1)+1:end);
% 畫 A_tdomain time plot
figure(7)
plot(h(look_mic, :), 'r');
hold on
plot(A_tdomain_forplot, 'b');
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('A\_tdomain')
legend('h', 'A\_tdomain')
xlabel('points')
ylabel('amplitude')
shg

error_A_tdomain = sum((A_tdomain_forplot - h(look_mic, 1:points_rir-hopsize*(osfac-1))).^2);

source_WPD = istft(S_WPD, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, ConjugateSymmetric=true, FrequencyRange='onesided');
source_WPD = source_WPD.';
source_max  = max(abs(source(1, :)));
source_WPD_max  = max(abs(source_WPD(1, :)));
ratio_source_WPD = source_max/source_WPD_max;
source_WPD = source_WPD.*ratio_source_WPD;
audiowrite('wav\source_WPD_all.wav', source_WPD(1, :), fs)

%%  dereverb with MINT (source_MINT) %%
A_tdomain_forMINT = A_tdomain(:, hopsize*(osfac-1)+1:end);
g_len = 6000;
weight_len = 1500;
dia_load_MINT = 10^(-7);
source_MINT = MINT(A_tdomain_forMINT, y, g_len, weight_len, dia_load_MINT);

% adjust source_predict 的最高點使之與 source 的一樣 %
source_max  = max(abs(source(1, :)));
source_MINT_max  = max(abs(source_MINT(1, :)));
ratio_source_MINT = source_max/source_MINT_max;
source_MINT = source_MINT.*ratio_source_MINT;

% 畫 source_predict time plot %
figure(7)
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
audiowrite('wav\source_all.wav', source(1, :), fs)

ratio_y_nodelay = 0.8 / max(abs(y_nodelay(look_mic, :))) ;
y_filemane_str = ['wav\y_nodelay_all-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_nodelay(look_mic, :)*ratio_y_nodelay, fs)

source_MINT_filemane_str = ['wav\source_predict_all_WPD_RLS_MINT-', string(reverberation_time), 'x', string(g_len), 'x', string(weight_len),...
    'x', string(dia_load_MINT) '.wav'];
source_MINT_filemane = join(source_MINT_filemane_str, '');
audiowrite(source_MINT_filemane, source_MINT(1, :), fs)

fprintf('done\n')
