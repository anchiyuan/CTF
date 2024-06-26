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
reverberation_time = 0.6;                                % Reverberation time (s)
points_rir = 12288;                                       % Number of rir points (需比 reverberation time 還長)
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 1;                                           % Disable high-pass filter

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

%% reconstruct RIR from H (h_reconstruct) %%
% h_reconstruct = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, H);
% 
% % 校正 h_reconstruct 之高度 %
% ratio_h_reconstruct = zeros(MicNum, 1);
% for i = 1:MicNum
%     ratio_h_reconstruct(i, :) = max(abs(h(i, :)))/max(abs(h_reconstruct(i, :)));
% end
% 
% h_reconstruct = h_reconstruct.*ratio_h_reconstruct;
% 
% h_reconstruct_forplot = h_reconstruct(:, hopsize*(osfac-1)+1:end);
% % 畫 reconstructed RIR time plot %
% figure(3)
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

%% 讀音檔 or 產生 white noise source (source) %%
Second = 30;
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

%% load y_wpe (Y_wpe) %%
y_wpe_filename_str = ['y_wpe_30-', string(reverberation_time), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);
y_wpe_transpose = y_wpe.';
[Y_wpe, ~, ~] = stft(y_wpe_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% calculate sigma %%
hight_view = 0.2;
sec_rough = ones(1, SorLen)*hight_view;
threshold_rough = 0.01;
for i = 1:SorLen
    if abs(source(:, i)) < threshold_rough
        sec_rough(:, i) = 0;
    end

end

sec_speech = ones(1, SorLen)*hight_view;
half_point_for_mean = 100;
threshold_speech = hight_view*half_point_for_mean/(half_point_for_mean*2+1);
for i = half_point_for_mean+1:SorLen-(half_point_for_mean+1)
    if mean(sec_rough(:, i-half_point_for_mean:i+half_point_for_mean)) < threshold_speech
        sec_speech(:, i) = 0;
    end

end

figure(20)
plot(source(1, :), 'r');
hold on
plot(sec_speech(1, :), 'b');
hold off
title('speech section')
xlabel('points')
ylabel('magnitude')
legend('source', 'section')
shg

forfac_Ryy = 0.999;
forfac_Ruu = 0.999;
Ryy = zeros(MicNum, MicNum, frequency);
Ruu = zeros(MicNum, MicNum, frequency);
for FrameNo = 1:NumOfFrame
    if mean(sec_speech(:, (FrameNo-1)*hopsize+1:(FrameNo-1)*hopsize+NFFT)) == 0
        for n = 1:frequency
            Ruu(:, :, n) = forfac_Ruu*Ruu(:, :, n) + (1-forfac_Ruu)*squeeze(Y_nodelay(n, FrameNo, :))*squeeze(Y_nodelay(n, FrameNo, :))';
        end

    else
        for n = 1:frequency
            Ryy(:, :, n) = forfac_Ryy*Ryy(:, :, n) + (1-forfac_Ryy)*squeeze(Y_nodelay(n, FrameNo, :))*squeeze(Y_nodelay(n, FrameNo, :))';
        end

    end

end

Rxx = Ryy - Ruu;
dia_load_W = 10^(0);
W = zeros(MicNum, MicNum, frequency);
X = zeros(frequency, NumOfFrame, MicNum);
for n = 1:frequency
    W(:, :, n) = inv(Ryy(:, :, n) + dia_load_W*eye(MicNum))*Rxx(:, :, n);
    for FrameNo = 1:NumOfFrame
        X(n, FrameNo, :) = W(:, :, n)'*squeeze(Y_nodelay(n, FrameNo, :));
    end

end

X = [X zeros(frequency, 1, MicNum)];

x = istft(X, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, ConjugateSymmetric=true, FrequencyRange='onesided');
x = x.';

% 畫 reconstructed RIR time plot %
figure(21)
plot(y_nodelay(look_mic, :), 'r');
hold on
plot(x(look_mic, :), 'b');
hold off
legend('y\_nodelay', 'x')
xlabel('points')
ylabel('amplitude')
shg

ratio_y_nodelay = 0.8 / max(abs(y_nodelay(look_mic, :))) ;
y_filemane_str = ['wav_WPD\y_nodelay_all-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_nodelay(look_mic, :)*ratio_y_nodelay, fs)

ratio_x = 0.8 / max(abs(x(look_mic, :))) ;
x_filemane_str = ['wav_WPD\x_all-', string(reverberation_time), '.wav'];
x_filemane = join(x_filemane_str, '');
audiowrite(x_filemane, x(look_mic, :)*ratio_x, fs)

%% WPD and Wiener (A S_WPD) %%
b = 2;
frequency_range = [800, 1500];
Lw = [6, 5, 3];    % 0–800, 800–1500, 1500–3000, and 3000–4000 Hz.
Lf = 1;
start_ini_frame = Lw(1) + L;
ini_frame = NumOfFrame;
reference_mic = 15;
diaload_R_bar = 1e3;
forfac_Rss = 0.999;
forfac_Rsy = 0.999;
dia_load_A = 10^(-7);
S_WPD = zeros(frequency, NumOfFrame);
A = zeros(MicNum, L, frequency);

% compute distance for a %
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((SorPos - MicPos(i, :)).^2));
end

% iteration process %
parfor n = 1 : frequency
    % choose Lw for different frequency %
    if (freqs_vector(n) <= frequency_range(1))
        Lw_choose = Lw(1);
    elseif (freqs_vector(n) <= frequency_range(2))
        Lw_choose = Lw(2);
    else
        Lw_choose = Lw(3);
    end
   
    % initialize Y_bar %
    Y_bar = zeros(MicNum*(Lw_choose-b+2), NumOfFrame)

    % initialize R_bar %
    R_bar = zeros(MicNum*(Lw_choose-b+2), MicNum*(Lw_choose-b+2));

    % construct a_bar %
    omega = 2*pi*freqs_vector(n);
    a = exp(-1j*omega/c*distance)./distance;
    a_bar = zeros(MicNum*(Lw_choose-b+2), 1);
    a_bar(1:MicNum, :) = a;

    for FrameNo = start_ini_frame:ini_frame
        % cascade right now processing frame Y_bar %
        Y_bar_now = squeeze(Y_nodelay(n, FrameNo, :));
        for frame = b:Lw_choose
            Y_bar_now = cat(1, Y_bar_now, squeeze(Y_nodelay(n, FrameNo-frame, :)));
        end

        Y_bar(:, FrameNo) = Y_bar_now;

        % compute sigma %
%         sigma = mean(abs(S(n, FrameNo)).^2, 'all');
%         sigma = mean(abs(Y_wpe(:, FrameNo, reference_mic)).^2, 'all');
        sigma = mean(abs(X(n, FrameNo-Lf:FrameNo+Lf, reference_mic)).^2, 'all');

        % accumulate R_bar %
        R_bar = R_bar + Y_bar_now*Y_bar_now'/sigma;
    end

    % compute WPD weight %
    w_bar = inv(R_bar + diaload_R_bar*eye(MicNum*(Lw_choose-b+2)))*a_bar/(a_bar'*inv(R_bar + diaload_R_bar*eye(MicNum*(Lw_choose-b+2)))*a_bar);
    S_WPD(n, :) = w_bar'*Y_bar;
    
    % compute Rss Rsy %
    S_WPD_frequencywise = w_bar'*Y_bar;
    Rss = zeros(L, L);
    Rsy = zeros(L, MicNum);
    for FrameNo= start_ini_frame+L-1:ini_frame
        S_choose = S_WPD_frequencywise(:, FrameNo-L+1:FrameNo).';
        Y_choose = squeeze(Y_delay(n, FrameNo, :));
        Rss = forfac_Rss*Rss + (1-forfac_Rss)*flip(S_choose) * flip(S_choose)';
        Rsy = forfac_Rsy*Rsy + (1-forfac_Rsy)*flip(S_choose) * Y_choose';
    end

    % compute A %
    A(:, :, n) = Rsy'/(Rss + dia_load_A*eye(L));
end

A_dB = mag2db(abs(A));
% 畫 A frequency plot %
figure(5);
mesh(L_vector, freqs_vector, squeeze(A_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A')
xlabel('frame')
ylabel('frequency(Hz)')
shg

% save source_WPD .wav 檔 %
source_WPD = istft(S_WPD, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, ConjugateSymmetric=true, FrequencyRange='onesided');
source_WPD = source_WPD.';
source_max  = max(abs(source(1, :)));
source_WPD_max  = max(abs(source_WPD(1, :)));
ratio_source_WPD = source_max/source_WPD_max;
source_WPD = source_WPD.*ratio_source_WPD;
source_WPD_filemane_str = ['wav_WPD\source_WPD_all-', string(reverberation_time), '.wav'];
source_WPD_filemane = join(source_WPD_filemane_str, '');
audiowrite(source_WPD_filemane, source_WPD(1, :), fs)
audiowrite('wav_WPD\source_all.wav', source(1, :), fs)

% S_WPD_dB = mag2db(abs(S_WPD));
% % 畫 S_WPD frequency plot %
% figure(6);
% mesh(1:1:size(S_WPD, 2), freqs_vector, S_WPD_dB)
% colorbar
% view(2)
% xlim([1 NumOfFrame])
% title('S\_WPD')
% xlabel('frame')
% ylabel('frequency(Hz)')
% shg

%% A 轉回時域 (A_tdomain) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
ratio_A_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_tdomain(i, :) = max(abs(h(i, :)))/max(abs(A_tdomain(i, hopsize*(osfac-1)+1:end)));
end

A_tdomain = A_tdomain.*ratio_A_tdomain;

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

% 算 matching error %
ATF = fft(h);
ATF_estimated = fft(A_tdomain(:, hopsize*(osfac-1)+1:end));
sum_norm = 0;
for i  = 1:MicNum
    norm_ATF = norm(ATF(i, :) - ATF_estimated(i, :));
    sum_norm = sum_norm + norm_ATF;
end

ME = sum_norm/MicNum;

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
% audiowrite('wav_WPD\source_all.wav', source(1, :), fs)
% 
ratio_y_nodelay = 0.8 / max(abs(y_nodelay(look_mic, :))) ;
y_filemane_str = ['wav_WPD\y_nodelay_all-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_nodelay(look_mic, :)*ratio_y_nodelay, fs)
% 
% source_MINT_filemane_str = ['wav_WPD\source_predict_all_WPD_Wiener_MINT-', string(reverberation_time), '.wav'];
% source_MINT_filemane = join(source_MINT_filemane_str, '');
% audiowrite(source_MINT_filemane, source_MINT(1, :), fs)

% save partial wav %
point_start_save = 25*fs;

audiowrite('wav_WPD\source_partial.wav', source(1, point_start_save:end), fs)

ratio_y_nodelay = 0.8 / max(abs(y_nodelay(look_mic, point_start_save:end))) ;
y_filemane_str = ['wav_WPD\y_nodelay_partial-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_nodelay(look_mic, point_start_save:end)*ratio_y_nodelay, fs)

ratio_y_wpe = 0.8 / max(abs(y_wpe(look_mic, point_start_save:end))) ;
y_filemane_str = ['wav_WPD\y_wpe_partial-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_wpe(look_mic, point_start_save:end)*ratio_y_nodelay, fs)

source_WPD_filemane_str = ['wav_WPD\source_WPD_partial-', string(reverberation_time), '.wav'];
source_WPD_filemane = join(source_WPD_filemane_str, '');
audiowrite(source_WPD_filemane, source_WPD(1, point_start_save:end), fs)

source_MINT_filemane_str = ['wav_WPD\source_predict_partial_WPD_Wiener_MINT-', string(reverberation_time), '.wav'];
source_MINT_filemane = join(source_MINT_filemane_str, '');
audiowrite(source_MINT_filemane, source_MINT(1, point_start_save:end), fs)

fprintf('done\n')
