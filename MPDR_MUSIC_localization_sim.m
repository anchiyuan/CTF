clc; clear;
close all;

tic

%% RIR parameter %%
SorNum = 1;                                              % source number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MicNum = 14;                                             % number of microphone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)

% sparse array %
MicPos = zeros(MicNum, 3);
sa_num = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
micnum_in_sa = 2;
sa_center = [1, 1, 1; 1, 4, 1; 4, 4, 1; 4, 1, 1];
sa_space = 0.01;
mic_space = 0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sparce_num = sa_num*micnum_in_sa;
for i = 1:sa_num
    sa_pos = zeros(micnum_in_sa, 3);
    for j = 1:micnum_in_sa/2
        sa_pos(j, :) = [sa_center(i, 1)-sa_space-(micnum_in_sa/2-j)*mic_space, sa_center(i, 2), sa_center(i, 3)];
        sa_pos(micnum_in_sa/2+j, :) = [sa_center(i, 1)+sa_space+(j-1)*mic_space, sa_center(i, 2), sa_center(i, 3) ];
    end

    MicPos(1+(i-1)*micnum_in_sa:i*micnum_in_sa, :) = sa_pos;
end

% MicPos(1:4, :) = [1, 1, 1; 1, 4, 1; 4, 4, 1; 4, 1, 1];

referencce_point = mean(sa_center);

% UCA %
UCA_num = MicNum - sparce_num;
radius_UCA = 0.25;
angle_UCA = 360/UCA_num;
for i = 1:UCA_num
    MicPos(i+sparce_num, :) = [referencce_point(:, 1)+radius_UCA*sind((i-1)*angle_UCA), referencce_point(:, 2)+radius_UCA*cosd((i-1)*angle_UCA), referencce_point(:, 3)];
end

SorPos = [1.5, 3.2, 1];                                    % source position (m)
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
angle_ground_truth2 = 360 - angle_ground_truth;
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

y_sparce = as(1:sparce_num, 1:SorLen);
y_UCA = as(sparce_num+1:end, 1:SorLen);
[Y_sparce, ~, ~] = stft(y_sparce.', fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
[Y_UCA, ~, ~] = stft(y_UCA.', fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
NumOfFrame = size(Y_sparce, 2);

%% compute Ryy %%
Ryy_sparce = zeros(sparce_num, sparce_num, frequency);
Ryy_UCA = zeros(UCA_num, UCA_num, frequency);
for n = 1:frequency
    for FrameNo = 1:NumOfFrame
        Ryy_sparce(:, :, n) = Ryy_sparce(:, :, n) + squeeze(Y_sparce(n, FrameNo, :))*squeeze(Y_sparce(n, FrameNo, :))';
        Ryy_UCA(:, :, n) = Ryy_UCA(:, :, n) + squeeze(Y_UCA(n, FrameNo, :))*squeeze(Y_UCA(n, FrameNo, :))';
    end
    
end

Ryy_sparce = Ryy_sparce/NumOfFrame;
Ryy_UCA = Ryy_UCA/NumOfFrame;



%% compute P_N %%
P_N = zeros(sparce_num, sparce_num, frequency);
for n = 1:frequency
   [U, D] = eig(Ryy_sparce(:, :, n));
   [~, D_max_index] = max(diag(D));
   U_S = U(:, D_max_index);
   P_N(:, :, n) = eye(sparce_num) - U_S*U_S';
end

%% redefine MicPos %%
MicPos = MicPos - referencce_point;

%% check angle function is convex or not %%
angle = 0:1:359;
array_output_power = zeros(frequency, size(angle, 2));
output_power_sum = zeros(size(angle, 2), 1);

dia_load_beamformer = 10^(-2);
frequency_lower_bound = 1;

angle_count = 0;

for ang = angle
    angle_count = angle_count + 1;
    kappa = [sind(ang), cosd(ang), 0];
    for n = frequency_lower_bound:frequency
        omega = 2*pi*freqs_vector(n);
        steer_vec = exp(1j*omega/c*kappa*MicPos(sparce_num+1:end, :).').';

        array_output_power(n, angle_count) = 1/(steer_vec'*inv(Ryy_UCA(:, :, n)+dia_load_beamformer*eye(UCA_num))*steer_vec);
        output_power_sum(angle_count, :) = output_power_sum(angle_count, :) + abs(array_output_power(n, angle_count));
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
    
    left_insert_output = GSS_MPDR_angle_obj_func(left_insert, frequency_lower_bound, frequency, freqs_vector, c, MicPos(sparce_num+1:end, :), Ryy_UCA, UCA_num, dia_load_beamformer);
    right_insert_output = GSS_MPDR_angle_obj_func(right_insert, frequency_lower_bound, frequency, freqs_vector, c, MicPos(sparce_num+1:end, :), Ryy_UCA, UCA_num, dia_load_beamformer);

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
distance = 0.1:0.1:2.1;
array_output_power = zeros(frequency, size(distance, 2));
output_power_sum = zeros(size(distance, 2), 1);
frequency_lower_bound = 1;

distance_count = 0;

for dis = distance
    distance_count = distance_count + 1;
    source_pos = [dis*sind(angle_final), dis*cosd(angle_final), 0];

    r = zeros(sparce_num, 1);
    for i = 1:sparce_num
        r(i, :) =  sqrt(sum((source_pos - MicPos(i, :)).^2));
    end

    for n = frequency_lower_bound:frequency
        omega = 2*pi*freqs_vector(n);
        steer_vec = exp(-1j*omega/c*r)./r;

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
localmax = islocalmax(output_power_sum);
output_power_sum(~localmax) = 0;
[~, max_index] = max(output_power_sum);
left_bound = distance(:, max_index-1);
right_bound = distance(:, max_index+1);

iteration_times_distnce = 0;
left_move = 1;
right_move = 1;

while 1
    iteration_times_distnce = iteration_times_distnce + 1;
    left_insert = right_bound - search_ratio*(right_bound-left_bound);
    right_insert = left_bound + search_ratio*(right_bound-left_bound);

    if right_move == 1 && left_move == 1
        left_insert_output = GSS_MUSIC_distance_obj_func(left_insert, angle_final, sparce_num, MicPos(1:sparce_num, :), frequency_lower_bound, frequency, freqs_vector, c, P_N);
        right_insert_output = GSS_MUSIC_distance_obj_func(right_insert, angle_final, sparce_num, MicPos(1:sparce_num, :), frequency_lower_bound, frequency, freqs_vector, c, P_N);
        right_move = 0;
        leftt_move = 0;

    elseif right_move == 1 && left_move == 0
        right_insert_output = left_insert_output;
        left_insert_output = GSS_MUSIC_distance_obj_func(left_insert, angle_final, sparce_num, MicPos(1:sparce_num, :), frequency_lower_bound, frequency, freqs_vector, c, P_N);
        right_move = 0;

    elseif right_move == 0 && left_move == 1
        left_insert_output = right_insert_output;
        right_insert_output = GSS_MUSIC_distance_obj_func(right_insert, angle_final, sparce_num, MicPos(1:sparce_num, :), frequency_lower_bound, frequency, freqs_vector, c, P_N);
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
