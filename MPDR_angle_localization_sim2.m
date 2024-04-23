clc; clear;
close all;

tic

%% RIR parameter %%
SorNum = 1;                                              % source number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MicNum = 20;                                             % number of microphone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)

% subarray %
sa_center = [1, 1, 1];
MicPos = zeros(MicNum, 3);
sa_space = 0.25;
mic_space = 0.02;
for i = 1:MicNum/2
    MicPos(i, :) = [sa_center(1, 1)-sa_space-(MicNum/2-i)*mic_space, sa_center(1, 2), sa_center(1, 3)];
    MicPos(MicNum/2+i, :) = [sa_center(1, 1)+sa_space+(i-1)*mic_space, sa_center(1, 2), sa_center(1, 3) ];
end

referencce_point = [2.5, 2.5, 1];

SorPos = [2.8, 2.9, 1];                                    % source position (m)
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
plot3(MicPos(:, 1), MicPos(:, 2), MicPos(:, 3), 'r.', 'MarkerSize', 15)
hold on
plot3(SorPos(:, 1), SorPos(:, 2), SorPos(:, 3), '*', 'MarkerSize', 20)
hold on
plot3(referencce_point(:, 1), referencce_point(:, 2), referencce_point(:, 3), 'x', 'MarkerSize', 15)
hold off
xlabel('x\_axis')
ylabel('y\_axis')
zlabel('z\_axis')
title('空間圖')
shg

angle_ground_truth = acos((SorPos-referencce_point)*[0; 1; 0]/norm(SorPos-referencce_point))/pi*180;
distance_ground_truth = norm(SorPos-referencce_point);

%% generate ground-truth RIR (h) %%
% 產生 RIR 和存.mat 檔 %
h = rir_generator(c, fs, MicPos, SorPos, room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
rir_filename_str = ['h\h_for_localization_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
save(rir_filemane, 'h')

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

[~, argmax_tf] = max(abs(h.'));

%% window parameter %%
NFFT = 1024;
hopsize = 256;
win = hamming(NFFT);
frequency = NFFT/2 + 1;
freqs_vector = linspace(0, fs/2, frequency);

%% 讀音檔 or 產生 white noise source (source) %%
Second = 23;
SorLen =  Second*fs;

% load speech source %
% [source_transpose, fs] = audioread('245.wav', [1, SorLen]);    % speech source
% source = source_transpose.';

% load white noise source %
source = wgn(1, SorLen, 0);                                    % white noise source

%% 產生麥克風訊號，先在時域上 convolution 再做 stft (y_nodelay y_delay Y_delay ) %%
% convolution source and RIR %
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(h(i, :), source);
end

y_nodelay = as(:, 1:SorLen);
y_nodelay_transpose = y_nodelay.';
[Y_nodelay, ~, ~] = stft(y_nodelay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
NumOfFrame = size(Y_nodelay, 2);

%% compute Ryy %%
Ryy = zeros(MicNum, MicNum, frequency);
for_fac = 0.001;
for n = 1:frequency
    for FrameNo = 1:NumOfFrame
%         Ryy(:, :, n) = for_fac*Ryy(:, :, n) + (1-for_fac)*squeeze(Y_nodelay(n, FrameNo, :))*squeeze(Y_nodelay(n, FrameNo, :))';
        Ryy(:, :, n) = Ryy(:, :, n) + squeeze(Y_nodelay(n, FrameNo, :))*squeeze(Y_nodelay(n, FrameNo, :))';

    end
    
end

Ryy = Ryy/NumOfFrame;

%% compute P_N %%
P_N = zeros(MicNum, MicNum, frequency);
for n = 1:frequency
   [U, D] = eig(Ryy(:, :, n));
   [~, D_max_index] = max(diag(D));
   U_S = U(:, D_max_index);
   P_N(:, :, n) = eye(MicNum) - U_S*U_S';

end

%% redefine MicPos %%
MicPos = MicPos - referencce_point;

%% check angle function is convex or not %%
angle = 0:1:359;
array_output_power = zeros(frequency, size(angle, 2));
output_power_sum = zeros(size(angle, 2), 1);

dia_load_beamformer = 10^(-2);
frequency_lower_bound = 400;

angle_count = 0;

for ang = angle
    angle_count = angle_count + 1;
    kappa = [sind(ang), cosd(ang), 0];
    for n = frequency_lower_bound:frequency
        omega = 2*pi*freqs_vector(n);
        steer_vec = exp(1j*omega/c*kappa*MicPos.').';

        array_output_power(n, angle_count) = 1/(steer_vec'*inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*steer_vec);
        output_power_sum(angle_count, :) = output_power_sum(angle_count, :) + abs(array_output_power(n, angle_count));

%         array_output_power(n, angle_count) = 1/(steer_vec'*P_N(:, :, n)*steer_vec);
%         output_power_sum(angle_count, :) = output_power_sum(angle_count, :) + abs(array_output_power(n, angle_count));

    end

end

[meshgrid_x, meshgrid_y] = meshgrid(angle ,freqs_vector(frequency_lower_bound:frequency));
figure(3);
mesh(meshgrid_x, meshgrid_y, 10*log10(abs(array_output_power(frequency_lower_bound:frequency, :))))
colorbar
view(2)
title('output power')
xlabel('angle')
ylabel('frequency')

figure(4)
plot(angle, output_power_sum.');
title('output power as function of angle')
xlabel('angle')
ylabel('output power')
shg

%% GSS for angle-wise freefield plane wave localization %%
[~, max_index] = max(output_power_sum);
left_bound = angle(:, max_index-1);
right_bound = angle(:, max_index+1);

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

angle_final = (left_bound+right_bound)/2;

%% check distance function is convex or not %%
distance = 0:0.1:2.1;
array_output_power = zeros(frequency, size(distance, 2));
output_power_sum = zeros(size(distance, 2), 1);
frequency_lower_bound = 1;

distance_count = 0;

for dis = distance
    distance_count = distance_count + 1;
    source_pos = [dis*sind(angle_ground_truth), dis*cosd(angle_ground_truth), 0];

    r = zeros(MicNum, 1);
    for i = 1 : MicNum
        r(i, :) =  sqrt(sum((source_pos - MicPos(i, :)).^2));
    end

    for n = frequency_lower_bound:frequency
        omega = 2*pi*freqs_vector(n);
        steer_vec = exp(-1j*omega/c*r)./r;

%         array_output_power(n, distance_count) = steer_vec'*Ryy(:, :, n)*steer_vec;
%         output_power_sum(distance_count, :) = output_power_sum(distance_count, :) + abs(array_output_power(n, distance_count));

        array_output_power(n, distance_count) = 1/(steer_vec'*P_N(:, :, n)*steer_vec);
        output_power_sum(distance_count, :) = output_power_sum(distance_count, :) + abs(array_output_power(n, distance_count));
        
    end

end

[meshgrid_x, meshgrid_y] = meshgrid(distance ,freqs_vector(frequency_lower_bound:frequency));
figure(5);
mesh(meshgrid_x, meshgrid_y, 10*log10(abs(array_output_power(frequency_lower_bound:frequency, :))))
colorbar
view(2)
title('pseudo spectrum')
xlabel('distance')
ylabel('frequency')


figure(6)
plot(distance, output_power_sum.');
title('spectrum magnitude as function of distance')
xlabel('distance')
ylabel('magnitude')
shg

%% GSS for distnce-wise MUSIC %%
output_power_sum = islocalmax(output_power_sum);
[~, max_index] = max(output_power_sum);
left_bound = distance(:, max_index-1);
right_bound = distance(:, max_index+1);

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
        left_insert_output = GSS_MUSIC_distance_obj_func(left_insert, angle_ground_truth, MicNum, MicPos, frequency_lower_bound, frequency, freqs_vector, c, P_N);
        right_insert_output = GSS_MUSIC_distance_obj_func(right_insert, angle_ground_truth, MicNum, MicPos, frequency_lower_bound, frequency, freqs_vector, c, P_N);
        right_move = 0;
        leftt_move = 0;

    elseif right_move == 1 && left_move == 0
        right_insert_output = left_insert_output;
        left_insert_output = GSS_MUSIC_distance_obj_func(left_insert, angle_ground_truth, MicNum, MicPos, frequency_lower_bound, frequency, freqs_vector, c, P_N);
        right_move = 0;

    elseif right_move == 0 && left_move == 1
        left_insert_output = right_insert_output;
        right_insert_output = GSS_MUSIC_distance_obj_func(right_insert, angle_ground_truth, MicNum, MicPos, frequency_lower_bound, frequency, freqs_vector, c, P_N);
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

distance_final = (left_bound+right_bound)/2;

toc
