clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

%% RIR parameter %%
SorNum = 2;                                              % source number
MicNum_TDOA = 8;                                         % TDOA麥克風數量
MicNum = 38;                                             % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)

% distributed 8 mic %
mic_x = [ 200 ; 300 ; 300 ; 200 ; 200 ; 300 ; 300 ; 200 ]./100;
mic_y = [ 200 ; 200 ; 300 ; 300 ; 200 ; 200 ; 300 ; 300 ]./100;
mic_z = [ 100 ; 100 ; 100 ; 100 ; 200 ; 200 ; 200 ; 200 ]./100;
MicPos = [mic_x, mic_y, mic_z,];

% ULA 30 mics %
MicStart = [210, 200, 100]/100;
spacing = 0.02;
for i = MicNum_TDOA+1:MicNum
    MicPos(i, :) = [MicStart(1, 1)+(i-(MicNum_TDOA+1))*spacing, MicStart(1, 2), MicStart(1, 3)];
end

SorPos = [210, 215, 110; 290, 290, 190]/100;                            % source position (m)
room_dim = [500, 600, 250]/100;                          % Room dimensions [x y z] (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reverberation_time = 0.6;                                % Reverberation time (s)
points_rir = 12288;                                       % Number of rir points (需比 reverberation time 還長)
look_mic = 10;
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
% % 產生 RIR 和存.mat 檔 %
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

% 畫 ground-truth RIR time plot %
figure(2)
subplot(2, 1, 1);
plot(squeeze(h(look_mic, 1, :)).', 'r');
title('source RIR')
xlabel('points')
ylabel('amplitude')

subplot(2, 1, 2);
plot(squeeze(h(look_mic, 2, :)).', 'r');
title('interferer RIR')
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

%% 讀音檔 (source) %%
Second = 23;
SorLen =  Second*fs;

% load speech source %
source_transpose = zeros(SorLen, SorNum);
[source_transpose(:, 1), ~] = audioread('245.wav', [1, SorLen]);
[source_transpose(:, 2), ~] = audioread('313.wav', [1, SorLen]);% speech source
source = source_transpose.';

%% compute source signal for frequency (S) %%
% source 轉頻域 %
source_transpose = source.';
[S, ~, ~] = stft(source_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

NumOfFrame = size(S, 2);
NumOfFrame_vector = 1:1:NumOfFrame;

%% 產生麥克風訊號，先在時域上 convolution 再做 stft %%
% convolution source and RIR %
hs_source = zeros(MicNum, points_rir+SorLen-1);
hs_interferer = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    hs_source(i, :) = conv(squeeze(h(i, 1, :)).', source(1, :));
    hs_interferer(i, :) = conv(squeeze(h(i, 2, :)).', source(2, :));
end

hs_noisy = hs_source + hs_interferer;

% 從 hs 擷取需要的麥克風信號 y %
extra_delay = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF

y_source_delay = zeros(MicNum, SorLen);
y_source_delay(:, extra_delay+1:end) = hs_source(:, 1:SorLen-extra_delay);
y_source_nodelay = hs_source(:, 1:SorLen);

y_interferer_delay = zeros(MicNum, SorLen);
y_interferer_delay(:, extra_delay+1:end) = hs_interferer(:, 1:SorLen-extra_delay);
y_interferer_nodelay = hs_interferer(:, 1:SorLen);

y_noisy = hs_noisy(:, 1:SorLen);

% 各種 y 轉頻域 %
[Y_source_delay, ~, ~] = stft(y_source_delay.', fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
[Y_interferer_delay, ~, ~] = stft(y_interferer_delay.', fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');


%% WPE %%
% % do wpe %
% y_source_wpe = wpe(y_source_nodelay.', 'wpe_parameter.m');
% y_source_wpe = y_source_wpe.';
% 
% y_interferer_wpe = wpe(y_interferer_nodelay.', 'wpe_parameter.m');
% y_interferer_wpe = y_interferer_wpe.';
% 
% % 存 wpe mat %
% y_wpe_filename_str = ['y_TIKR\y_source_wpe-', string(reverberation_time), '.mat'];
% y_wpe_filename = join(y_wpe_filename_str, '');
% save(y_wpe_filename, 'y_source_wpe')
% 
% y_wpe_filename_str = ['y_TIKR\y_interferer_wpe-', string(reverberation_time), '.mat'];
% y_wpe_filename = join(y_wpe_filename_str, '');
% save(y_wpe_filename, 'y_interferer_wpe')

% load y_wpe %
y_wpe_filename_str = ['y_TIKR\y_source_wpe-', string(reverberation_time), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

y_wpe_filename_str = ['y_TIKR\y_interferer_wpe-', string(reverberation_time), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

%% DAS beamformer %%
% 算 mic 與 source 之距離 %
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, 1) =  sqrt(sum((SorPos(1, :) - MicPos(i, :)).^2));
    distance(i, 2) =  sqrt(sum((SorPos(2, :) - MicPos(i, :)).^2));
end

% 算 freefield ATF %
a = zeros(MicNum, SorNum, frequency);
for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, 1, n) = exp(-1j*omega/c*distance(:, 1))./distance(:, 1);
    a(:, 2, n) = exp(-1j*omega/c*distance(:, 2))./distance(:, 2);
end

% 算 DAS weight %
w = a/MicNum;

% 算 Y_DAS %
[Y_source_wpe, ~, ~] = stft(y_source_wpe.', fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
[Y_interferer_wpe, ~, ~] = stft(y_interferer_wpe.', fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

S_source_DAS = zeros(frequency, NumOfFrame);
S_interferer_DAS = zeros(frequency, NumOfFrame);
for FrameNo= 1:NumOfFrame
    for n = 1: frequency
         S_source_DAS(n, FrameNo) = w(:, 1, n)'*squeeze(Y_source_wpe(n, FrameNo, :));
         S_interferer_DAS(n, FrameNo) = w(:, 2, n)'*squeeze(Y_interferer_wpe(n, FrameNo, :));
    end  

end

%% predict CTF with Kalman filter (A) %%
start_ini_frame = L;
ini_frame = NumOfFrame;

% Kalman filter %
A_source = zeros(MicNum, L, frequency);
A_interferer = zeros(MicNum, L, frequency);
tic
parfor i = 1:MicNum
    for n =1:frequency
        weight_source = zeros(L, 1);
        weight_interferer = zeros(L, 1);
        P_source = 0.5*eye(L);    % error covariance matrix %
        P_interferer = 0.5*eye(L);    % error covariance matrix %
        K_source = zeros(L, 1);    % Kalman gain %
        K_interferer = zeros(L, 1);    % Kalman gain %
        R_source = 10^(-3);    % measurement noise covariance matrix %
        R_interferer = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_ini_frame:ini_frame
            % no time update only have measurement update %
            K_source = P_source*flip(S_source_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(S_source_DAS(n, FrameNo-L+1:FrameNo)))*P_source*flip(S_source_DAS(n, FrameNo-L+1:FrameNo).') + R_source);
            K_interferer = P_interferer*flip(S_interferer_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(S_interferer_DAS(n, FrameNo-L+1:FrameNo)))*P_interferer*flip(S_interferer_DAS(n, FrameNo-L+1:FrameNo).') + R_interferer);
            weight_source = weight_source + K_source*(conj(Y_source_delay(n, FrameNo, i)) - conj(flip(S_source_DAS(n, FrameNo-L+1:FrameNo)))*weight_source);
            weight_interferer = weight_interferer + K_interferer*(conj(Y_interferer_delay(n, FrameNo, i)) - conj(flip(S_interferer_DAS(n, FrameNo-L+1:FrameNo)))*weight_interferer);
            P_source = P_source - K_source*conj(flip(S_source_DAS(n, FrameNo-L+1:FrameNo)))*P_source;
            P_interferer = P_interferer - K_interferer*conj(flip(S_interferer_DAS(n, FrameNo-L+1:FrameNo)))*P_interferer;
        end
    
        A_source(i, :, n) = weight_source';
        A_interferer(i, :, n) = weight_interferer';
    end

end
toc

%% A 轉回時域且算 NRMSPM (A_tdomain NRMSPM) %%
A_source_forrecon = zeros(frequency, L, MicNum);
A_interferer_forrecon = zeros(frequency, L, MicNum);

for i = 1 : MicNum
    A_source_forrecon(:, :, i) = squeeze(A_source(i, :, :)).';
    A_interferer_forrecon(:, :, i) = squeeze(A_interferer(i, :, :)).';
end

A_source_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_source_forrecon);
A_interferer_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_interferer_forrecon);
A_source_tdomain = A_source_tdomain(:, hopsize*(osfac-1)+1:end);
A_interferer_tdomain = A_interferer_tdomain(:, hopsize*(osfac-1)+1:end);
ratio_A_source_tdomain = zeros(MicNum, 1);
ratio_A_interferer_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_source_tdomain(i, :) = max(abs(h(i, 1, :)))/max(abs(A_source_tdomain(i, :)));
    ratio_A_interferer_tdomain(i, :) = max(abs(h(i, 2, :)))/max(abs(A_interferer_tdomain(i, :)));
end

A_source_tdomain = A_source_tdomain.*ratio_A_source_tdomain;
A_interferer_tdomain = A_interferer_tdomain.*ratio_A_interferer_tdomain;

% 畫 A_tdomain time plot
figure(3)
subplot(2, 1, 1);
plot(squeeze(h(look_mic, 1, :)).', 'r');
hold on
plot(A_source_tdomain(look_mic, :), 'b');
hold off
title('estimated source RIR')
legend('ground-truth RIR', 'estimated RIR')
xlabel('points')
ylabel('amplitude')
shg

subplot(2, 1, 2);
plot(squeeze(h(look_mic, 2, :)).', 'r');
hold on
plot(A_interferer_tdomain(look_mic, :), 'b');
hold off
title('estimated interferer RIR')
legend('ground-truth RIR', 'estimated RIR')
xlabel('points')
ylabel('amplitude')
shg

%  算 NRMSPM %
h_NRMSPM = reshape(squeeze(h(:, 1, :)).', [MicNum*points_rir 1]);
aa_NRMSPM = reshape(A_source_tdomain.', [MicNum*points_rir 1]);
NRMSPM_source = 20*log10(norm(h_NRMSPM-h_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(h_NRMSPM));

h_NRMSPM = reshape(squeeze(h(:, 2, :)).', [MicNum*points_rir 1]);
aa_NRMSPM = reshape(A_interferer_tdomain.', [MicNum*points_rir 1]);
NRMSPM_interferer = 20*log10(norm(h_NRMSPM-h_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(h_NRMSPM));

%%  speech enhancement with MPDR (source_MPDR) %%
MPDR_mode = 'freefield_ATF';    % 'ATF' 'RTF' 'freefield_ATF'
if strcmp(MPDR_mode, 'ATF')
    % generate estimated ATF %
    frequency_ATF = points_rir/2 + 1;
    ATF = fft(A_source_tdomain, points_rir, 2);
    ATF = ATF(:, 1:frequency_ATF);

elseif strcmp(MPDR_mode, 'RTF')
    % generate estimated RTF %
    frequency_ATF = points_rir/2 + 1;
    ATF = fft(A_source_tdomain, points_rir, 2);
    ATF = ATF(:, 1:frequency_ATF);
    ATF = ATF./ATF(1, :);

elseif strcmp(MPDR_mode, 'freefield_ATF')
    % generate freefield ATF %
    frequency_ATF = points_rir/2 + 1;
    ATF = zeros(MicNum, frequency_ATF);
    frequency_ATF_vector = linspace(0, fs/2, frequency_ATF);
    for n = 1:frequency_ATF
        omega = 2*pi*frequency_ATF_vector(n);
        ATF(:, n) = exp(-1j*omega/c*distance(:, 1))./distance(:, 1);
    end

end

% noisy 麥克風訊號轉 stft %
[Y_noisy, ~, ~] = stft(y_noisy.', fs, Window=hamming(points_rir), OverlapLength=points_rir-points_rir/4, FFTLength=points_rir, FrequencyRange='onesided');
NumOfFrame_ATF =  size(Y_noisy, 2);

% compute Ryy %
Ryy = zeros(MicNum, MicNum, frequency_ATF);
for n = 1:frequency_ATF
    for FrameNo = 1:NumOfFrame_ATF
        Ryy(:, :, n) = Ryy(:, :, n) + squeeze(Y_noisy(n, FrameNo, :))*squeeze(Y_noisy(n, FrameNo, :))';
    end
    
end

Ryy = Ryy/NumOfFrame_ATF;

% compute MPDR weight %
dia_load_beamformer = 10^(2);    % critical for PESQ or SDR, its value depends on checking magnitude of Ryy
w_MPDR = zeros(MicNum, frequency);
for n = 1:frequency_ATF
    w_MPDR(:, n) = inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*ATF(:, n)/(ATF(:, n)'*inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*ATF(:, n));
end

% do MPDR %
S_MPDR = zeros(frequency_ATF, NumOfFrame_ATF);
for n = 1:frequency_ATF
    for FrameNo = 1:NumOfFrame_ATF
        S_MPDR(n, FrameNo) = w_MPDR(:, n)'*squeeze(Y_noisy(n, FrameNo, :));
    end
    
end

% ifft get source_MPDR %
source_MPDR_transpose = istft(S_MPDR, fs, Window=hamming(points_rir), OverlapLength=points_rir-points_rir/4, FFTLength=points_rir, ConjugateSymmetric=true, FrequencyRange='onesided');
source_MPDR = source_MPDR_transpose.';

%% save .wav 檔 %%
point_start_save = 18*fs;

audiowrite('wav_TIKR\source_ground-truth.wav', source(1, point_start_save:end), fs)
audiowrite('wav_TIKR\interferer_ground-truth.wav', source(2, point_start_save:end), fs)

ratio_y_noisy = 0.8 / max(abs(y_noisy(look_mic, point_start_save:end)));
y_filemane_str = ['wav_TIKR\y_noisy_', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y_noisy(look_mic, point_start_save:end)*ratio_y_noisy, fs)

ratio_source_MPDR = 1 / max(abs(source_MPDR(:, point_start_save:end)));
source_MPDR_filemane_str = ['wav_TIKR\source_MPDR_', string(MPDR_mode), '_', string(reverberation_time), '.wav'];
source_MPDR_filemane = join(source_MPDR_filemane_str, '');
audiowrite(source_MPDR_filemane, source_MPDR(1, point_start_save:end)*ratio_source_MPDR, fs)
