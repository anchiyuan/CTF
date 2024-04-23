clc; clear;
close all;

tic

%% parameters setting %%
SorNum = 1;
c = 343;
Fs = 48000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 16000;           % 欲 resample 成的取樣頻率
MicNum = 6;           % 實驗麥克風數量
points_rir = 2048;    % 自行設定想要輸出的 RIR 長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ULA %
MicStart = [0, 0, 0];
spacing = 0.07;
MicPos = zeros(MicNum, 3);
for i = 1:MicNum
    MicPos(i, :) = [MicStart(1, 1) + (i-1)*spacing, MicStart(1, 2), MicStart(1, 3)];
end

%% stft parameter %%
NFFT = 1024;
hopsize = 256;

% windows %
win = hamming(NFFT);
osfac = round(NFFT/hopsize);

frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);    % (len(win) + len(win) - 1) + points_rir - 1
L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);

%% read mic 音檔再做 stft (Y) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Second = 28;             % 使用時間長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*Fs;
sorLen =  Second*fs;

% load y %
y = zeros(MicNum, SorLen);
for i = 1:MicNum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_filename_str = ['wav_exp\', string(i), '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_filename = join(y_filename_str, '');
    [y(i, :), ~] = audioread(y_filename, [1, SorLen]);
end

% resample %
y = resample(y, 1, Fs/fs, Dimension=2);

% y do stft get Y %
y_transpose = y.';
[Y, ~, ~] = stft(y_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
NumOfFrame = size(Y, 2);

% delay y_nodelay to get y_delay %
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, sorLen);
y_delay(:, extra_delay_y+1:end) = y(:, 1:sorLen-extra_delay_y);

% y_delay 轉頻域 to get Y_delay %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPE (Y_wpe) %%
% % do wpe %
% y_wpe = wpe(y.', 'wpe_parameter.m');
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

% y_wpe 轉頻域 %
y_wpe_transpose = y_wpe.';
[Y_wpe, ~, ~] = stft(y_wpe_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% compute Ryy %%
Ryy = zeros(MicNum, MicNum, frequency);
for n = 1:frequency
    for FrameNo = 1:NumOfFrame
        Ryy(:, :, n) = Ryy(:, :, n) + squeeze(Y(n, FrameNo, :))*squeeze(Y(n, FrameNo, :))';
    end
    
end

Ryy = Ryy/NumOfFrame;

%% check angle function is convex or not %%
angle = -90:1:90;
output_abs_angle = zeros(size(angle, 2), 1);
dia_load_beamformer = 10^(-2);
angle_count = 0;

for ang = angle
    angle_count = angle_count + 1;
    kappa = [sind(ang), cosd(ang), 0];
    for n = 1:frequency
        omega = 2*pi*freqs_vector(n);
        steer_vec = exp(1j*omega/c*kappa*MicPos.').';
        array_output_power = 1/(steer_vec'*inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*steer_vec);
        output_abs_angle(angle_count, :) = output_abs_angle(angle_count, :) + abs(array_output_power);

    end

end

figure(1)
plot(angle, output_abs_angle.');
title('angle function')
xlabel('angle')
ylabel('output power')
shg

%% GSS for angle-wise freefield plane wave localization %%
[~, max_integer_angle] = max(output_abs_angle);
left_bound = angle(:, max_integer_angle-1);
right_bound = angle(:, max_integer_angle+1);
golden_ratio = (1+sqrt(5))/2;
search_ratio = golden_ratio - 1;
stop_criteron = 1e-7;

iteration_times_angle = 0;
while 1
    iteration_times_angle = iteration_times_angle + 1;
    left_insert = right_bound - search_ratio*(right_bound-left_bound);
    right_insert = left_bound + search_ratio*(right_bound-left_bound);
    
    left_insert_output = GSS_MPDR_angle_obj_func(left_insert, frequency, freqs_vector, c, MicPos, Ryy, MicNum, dia_load_beamformer);
    right_insert_output = GSS_MPDR_angle_obj_func(right_insert, frequency, freqs_vector, c, MicPos, Ryy, MicNum, dia_load_beamformer);

    if left_insert_output > right_insert_output
        right_bound = right_insert;
    elseif left_insert_output <= right_insert_output
        left_bound = left_insert;
    end

    if right_bound-left_bound < stop_criteron
        break
    end

end

final_angle = (left_bound+right_bound)/2;

%% GSS for distnce-wise RTF cosine similarity %%
mode = 'spherical wave';    % 'plane wave' or 'spherical wave'
left_bound = 0;
right_bound = 2;
golden_ratio = (1+sqrt(5))/2;
search_ratio = golden_ratio - 1;
stop_criteron = 1e-7;

iteration_times_distnce = 0;
left_move = 1;
right_move = 1;
while 1
    iteration_times_distnce = iteration_times_distnce + 1;
    left_insert = right_bound - search_ratio*(right_bound-left_bound);
    right_insert = left_bound + search_ratio*(right_bound-left_bound);

    if right_move == 1 && left_move == 1
        left_insert_output = GSS_CTF_distance_obj_func(left_insert, final_angle, MicNum, SorNum, MicPos, fs, c, NFFT, hopsize, points_rir, Y_wpe, Y_delay, mode);
        right_insert_output = GSS_CTF_distance_obj_func(right_insert, final_angle, MicNum, SorNum, MicPos, fs, c, NFFT, hopsize, points_rir, Y_wpe, Y_delay, mode);
        right_move = 0;
        leftt_move = 0;

    elseif right_move == 1 && left_move == 0
        right_insert_output = left_insert_output;
        left_insert_output = GSS_CTF_distance_obj_func(left_insert, final_angle, MicNum, SorNum, MicPos, fs, c, NFFT, hopsize, points_rir, Y_wpe, Y_delay, mode);
        right_move = 0;

    elseif right_move == 0 && left_move == 1
        left_insert_output = right_insert_output;
        right_insert_output = GSS_CTF_distance_obj_func(right_insert, final_angle, MicNum, SorNum, MicPos, fs, c, NFFT, hopsize, points_rir, Y_wpe, Y_delay, mode);
        leftt_move = 0;

    end

    if left_insert_output > right_insert_output
        right_bound = right_insert;
        right_move = 1;

    elseif left_insert_output <= right_insert_output
        left_bound = left_insert;
        left_move = 1;

    end

    if right_bound-left_bound < stop_criteron
        break

    end

end

final_distance = (left_bound+right_bound)/2;

toc
