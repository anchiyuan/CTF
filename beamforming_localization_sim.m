clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

tic

%% RIR parameter %%
SorNum = 1;                                              % source number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MicNum = 6;                                             % number of microphone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
points_rir = 2048;                                       % Number of rir points (需比 reverberation time 還長)
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

referencce_point = [1.05, 1.5, 1];
angle_ground_truth = acos((SorPos-referencce_point)*([0; 1; 0])/norm(SorPos-referencce_point))/pi*180;
distance_ground_truth = norm(SorPos-referencce_point);

%% generate ground-truth RIR (h) %%
% % 產生 RIR 和存.mat 檔 %
% h = rir_generator(c, fs, MicPos, SorPos, room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
% rir_filename_str = ['h\h_for_localization_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
% rir_filemane = join(rir_filename_str, '');
% save(rir_filemane, 'h')

% load RIR 的 .mat 檔 %
rir_filename_str = ['h\h_for_localization_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

look_mic = 1;
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
% 畫 ground-truth RIR time plot %
figure(2)
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

% 轉頻域 %
y_nodelay_transpose = y_nodelay.';
[Y_nodelay, ~, ~] = stft(y_nodelay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPE (y_wpe) %%
% % do wpe %
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
% y_wpe = y_wpe.';
% 
% % 存 wpe mat %
% y_wpe_filename_str = ['y\y_wpe_for_localization-', string(reverberation_time), 'x', string(MicNum), '.mat'];
% y_wpe_filename = join(y_wpe_filename_str, '');
% save(y_wpe_filename, 'y_wpe')

% load y_wpe %
y_wpe_filename_str = ['y\y_wpe_for_localization-', string(reverberation_time), 'x', string(MicNum), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

% y_wpe 轉頻域 %
y_wpe_transpose = y_wpe.';
[Y_wpe, ~, ~] = stft(y_wpe_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% compute Ryy %%
Ryy = zeros(MicNum, MicNum, frequency);
for n = 1:frequency
    for FrameNo = 1:NumOfFrame
        Ryy(:, :, n) = Ryy(:, :, n) + squeeze(Y_nodelay(n, FrameNo, :))*squeeze(Y_nodelay(n, FrameNo, :))';
    end
    
end

Ryy = Ryy/NumOfFrame;

%% redefine MicPos %%
% ULA %
MicStart = [0-(MicNum-1)/2*spacing, 0, 0];
MicPos = zeros(MicNum, 3);
for i = 1:MicNum
    MicPos(i, :) = [MicStart(1, 1)+(i-1)*spacing, MicStart(1, 2), MicStart(1, 3)];
end

%% check angle function is convex or not %%
angle = -90:1:90;
output_abs_angle = zeros(size(angle, 2), 1);
dia_load_beamformer = 10^(-2);
angle_count = 0;
frequency_lower_bound = 400;

for ang = angle
    angle_count = angle_count + 1;
    kappa = [sind(ang), cosd(ang), 0];
    for n = frequency_lower_bound:frequency
        omega = 2*pi*freqs_vector(n);
        steer_vec = exp(1j*omega/c*kappa*MicPos.').';
        array_output_power = 1/(steer_vec'*inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*steer_vec);
        output_abs_angle(angle_count, :) = output_abs_angle(angle_count, :) + abs(array_output_power);

    end

end

figure(3)
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
    
    left_insert_output = GSS_MPDR_angle_obj_func(left_insert, frequency_lower_bound, frequency, freqs_vector, c, MicPos, Ryy, MicNum, dia_load_beamformer);
    right_insert_output = GSS_MPDR_angle_obj_func(right_insert, frequency_lower_bound, frequency, freqs_vector, c, MicPos, Ryy, MicNum, dia_load_beamformer);

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
mode = 'plane wave';    % 'plane wave' or 'spherical wave'
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
