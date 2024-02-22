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

% 畫 ground-truth RIR time plot %
figure(1)
look_mic = 10;
plot(h(look_mic, :), 'r');
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h')
xlabel('points')
ylabel('amplitude')
shg

%% compute ground-truth CTF (H) %%
NFFT = 1024;
hopsize = 256;

% windows %
osfac = round(NFFT/hopsize);
win = hamming(NFFT);
coe = zeros(hopsize, 1);
for i = 1:osfac
    coe = coe + win((i-1)*hopsize+1:i*hopsize).^2;
end
coe = repmat(coe, [osfac, 1]);
awin = win./sqrt(NFFT*coe);    % analysis window
swin = win./sqrt(coe/NFFT);    % synthesis window

% compute H %
frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);
H = zeros(frequency, L, MicNum);
for k = 1:frequency
    zeta = xcorr(awin,swin).*exp(1i*2*pi*(k-1)*(-NFFT+1:NFFT-1)'/NFFT);
    for i = 1:MicNum
        H_temp = conv(h(i, :), zeta.');
        H(k, :, i) = H_temp(:, hopsize:hopsize:end);
    end
end

L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);
H_dB = mag2db(abs(H));

% 畫 ground-truth RIR frequency plot %
figure(2);
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

% source = wgn(1, SorLen, 0);    % white noise source

%% compute source signal for frequency (S) %%
% source 轉頻域 %
S = my_STFT(source, NFFT, hopsize, win);
NumOfFrame = size(S, 2);
NumOfFrame_vector = 1:1:NumOfFrame;
freqs_vector = linspace(0, fs/2, frequency);
S_dB = mag2db(abs(S));

% 畫 source frequency plot %
figure(3);
mesh(NumOfFrame_vector, freqs_vector, S_dB)
colorbar
view(2)
xlim([1 NumOfFrame])
title('S')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% RIR mix source 先在時域上 convolution 再做 stft (y and Y) %%
% convolution source and RIR %c
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(h(i, :), source);
end

extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y = zeros(MicNum, SorLen);
y(:, extra_delay_y+1:end) = as(:, 1:SorLen-extra_delay_y);

% y 轉頻域 %
Y = zeros(frequency, NumOfFrame, MicNum);
for i = 1:MicNum
    Y(:, :, i) = my_STFT(y(i, :), NFFT, hopsize, win);
end

%% WPE %%
% y_nodelay = as(:, 1:SorLen);
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
% y_wpe = y_wpe.';
% 
% y_wpe_filename_str = ['y_all_wpe-', string(reverberation_time), '.mat'];
% y_wpe_filename = join(y_wpe_filename_str, '');
% save(y_wpe_filename, 'y_wpe')

y_wpe_filename_str = ['y_wpe-', string(reverberation_time), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

% ratio_y = 1 / max(abs(y_nodelay(look_mic, :))) ;
% ratio_y_wpe = 1 / max(abs(y_wpe(look_mic, :))) ;
% 
% figure(4);
% plot(y_nodelay(look_mic, :), 'r');
% hold on
% plot(y_wpe(look_mic, :), 'b');
% hold off
% ylim([-1.1 1.1])
% title('y and y\_wpe')
% xlabel('points')
% ylabel('magnitude')
% legend('y', 'y\_wpe')
% shg

%% DAS beamformer %%
a = zeros(MicNum, SorNum, frequency);
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((SorPos - MicPos(i, :)).^2));
end

for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, :, n) = exp(-1j*omega/c*distance)./distance;
end

w = a/MicNum;

Y_wpe = zeros(frequency, NumOfFrame, MicNum);
for i = 1:MicNum
    Y_wpe(:, :, i) = my_STFT(y_wpe(i, :), NFFT, hopsize, win);
end

Y_DAS = zeros(frequency, NumOfFrame);
for FrameNo = 1:NumOfFrame
    for n = 1: frequency
         Y_DAS(n, FrameNo) = w(:, :, n)'*squeeze(Y_wpe(n, FrameNo, :));
    end  

end

%% predict CTF with LMS or RLS (A_ini) %%
start_ini_frame = ceil((12*fs - NFFT)/hopsize) + 1;    % wpe 第12秒開始穩定
ini_frame = NumOfFrame;    % 一個 frame 約等於0.016秒  看要娶幾個來 initialize

% 矯正 Y_DAS 的 PSD 使之與 S 一樣 %
PSD_S = sum(abs(S(:, start_ini_frame:ini_frame)).^2, "all");
PSD_Y_DAS = sum(abs(Y_DAS(:, start_ini_frame:ini_frame)).^2, "all");
PSD_ratio = PSD_S/PSD_Y_DAS;
Y_DAS = Y_DAS*sqrt(PSD_ratio);

% LMS %
% A_ini= zeros(MicNum, L, frequency);
% error = zeros(MicNum, (ini_frame-start_ini_frame+1), frequency);
% learn_rate = 1e-7;
% for i = 1:MicNum
%     for n = 1: frequency
%         weight = zeros(L, 1);
%         for FrameNo = start_ini_frame:ini_frame
%             Y_DAS_choose = flip(Y_DAS(n, FrameNo-L+1:FrameNo).');
%             Y_choose = Y(n, FrameNo, i);
%             error(i, FrameNo-start_ini_frame+1, n) = Y_choose - weight'*Y_DAS_choose;
%             weight = weight + learn_rate*Y_DAS_choose*conj(error(i, FrameNo-start_ini_frame+1, n));
%         end
% 
%         A_ini(i, :, n) = weight';
%     end
% 
% end

% RLS %
A_ini = zeros(MicNum, L, frequency);
error = zeros(MicNum, (ini_frame-start_ini_frame+1), frequency);
lambda = 0.99;
for i = 1:MicNum
    for n = 1: frequency
        weight = zeros(L, 1);
        P = (1e-2)^-1*eye(L);
        for FrameNo = start_ini_frame:ini_frame
            Y_DAS_choose = flip(Y_DAS(n, FrameNo-L+1:FrameNo).');
            Y_choose = Y(n, FrameNo, i);
            error(i, FrameNo-start_ini_frame+1, n) = Y_choose - weight'*Y_DAS_choose;
            k_RLS = lambda^-1*P*Y_DAS_choose/(1 + lambda^-1*Y_DAS_choose'*P*Y_DAS_choose);
            weight = weight + k_RLS*conj(error(i, FrameNo-start_ini_frame+1, n));
            P = lambda^-1*P - lambda^-1*k_RLS*Y_DAS_choose'*P;
        end

        A_ini(i, :, n) = weight';
    end

end

error_abs_mean = mean(mean(abs(error), 1), 3);

% 畫 error curve %
figure(5)
plot(error_abs_mean)
title('mean absolute error')
xlabel('iteration')
ylabel('magnitude')

A_ini_dB = mag2db(abs(A_ini));

% 畫 A_ini frequency plot %
figure(6);
mesh(L_vector, freqs_vector, squeeze(A_ini_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A\_ini')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% A_ini 轉回時域 (A_ini_tdomain) %%
A_ini_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_ini_forplot(:, :, i) = squeeze(A_ini(i, :, :)).';
end

A_ini_tdomain = reconstruct_RIR(points_rir, NFFT, hopsize, L, win, frequency, MicNum, A_ini_forplot);
A_ini_tdomain_forplot = A_ini_tdomain(look_mic, hopsize*(osfac-1)+1:end);

% max_h = max(abs(h(look_mic, :)));
% max_A_ini_tdomain_forplot = max(abs(A_ini_tdomain_forplot));
% ratio_A_ini_tdomain_forplot = max_h/max_A_ini_tdomain_forplot;
% A_ini_tdomain_forplot = A_ini_tdomain_forplot*ratio_A_ini_tdomain_forplot;

figure(7)
plot(h(look_mic, :), 'r');
hold on
plot(A_ini_tdomain_forplot, 'b');
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('A\_ini\_tdomain')
legend('h', 'A\_ini\_tdomain')
xlabel('points')
ylabel('amplitude')
shg

error_A_ini_tdomain = sum((A_ini_tdomain_forplot - h(look_mic, 1:points_rir-hopsize*(osfac-1))).^2);

%%  dereverb with MINT (source_MINT) %%
A_ini_tdomain_forMINT = A_ini_tdomain(:, hopsize*(osfac-1)+1:end);
g_len = 6000;
weight_len = 1500;
dia_load_MINT = 10^(-7);
source_MINT = MINT(A_ini_tdomain_forMINT, y, g_len, weight_len, dia_load_MINT);

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
audiowrite('wav\source_all.wav', source(1, :), fs)

ratio_y = 0.8 / max(abs(y(look_mic, :))) ;
y_filemane_str = ['wav\y_all-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y(look_mic, :)*ratio_y, fs)

source_MINT_filemane_str = ['wav\source_predict_all_RLS_MINT-', string(reverberation_time), 'x', string(g_len), 'x', string(weight_len),...
    'x', string(dia_load_MINT) '.wav'];
source_MINT_filemane = join(source_MINT_filemane_str, '');
audiowrite(source_MINT_filemane, source_MINT(1, :), fs)

fprintf('done\n')
