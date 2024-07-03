clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

tic

%% RIR parameter %%
SorNum = 1;                                              % source number
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

SorPos = [210, 215, 110]/100;                            % source position (m)
room_dim = [500, 600, 250]/100;                          % Room dimensions [x y z] (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reverberation_time = 0.2;                                % Reverberation time (s)
points_rir = 4096;                                       % Number of rir points (需比 reverberation time 還長)
look_mic = 38;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 1;                                           % Disable high-pass filter

figure(1);
plot3( [0 room_dim(1, 1) room_dim(1, 1) 0 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1)], ...
       [0 0 room_dim(1, 2) room_dim(1, 2) 0 0 0 room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) 0 0 0], ...
       [0 0 0 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0] , 'k')
hold on
plot3(MicPos(:, 1), MicPos(:, 2), MicPos(:, 3), 'r.', 'MarkerSize', 15)
hold on
plot3(SorPos(:, 1), SorPos(:, 2), SorPos(:, 3), '*', 'MarkerSize', 20)
hold off
xlabel('x\_axis')
ylabel('y\_axis')
zlabel('z\_axis')
title('空間圖')
shg

%% load ground-truth RIR (h) %%
% % 產生 RIR 和存.mat 檔 %
% h = rir_generator(c, fs, MicPos, SorPos, room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
% rir_filename_str = ['h\h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
% rir_filemane = join(rir_filename_str, '');
% save(rir_filemane, 'h')

rir_filename_str = ['h\h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

%% window parameter %%
NFFT = 1024;
hopsize = 256;

win = hamming(NFFT);
osfac = round(NFFT/hopsize);

frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);    % (len(win) + len(win) - 1) + points_rir - 1
freqs_vector = linspace(0, fs/2, frequency);

%% 讀音檔 (source) %%
Second = 23;
SorLen =  Second*fs;

% load speech source %
[source_transpose, fs] = audioread('245.wav', [1, SorLen]);
source = source_transpose.';

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
NumOfFrame = size(Y_delay, 2);

%% load y_wpe (y_wpe) %%
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

%% ASPSO %%
start_ini_frame = L;
ini_frame = NumOfFrame;
A = zeros(MicNum, L, frequency);

tic
parfor n =1:frequency    % 每個頻率獨立
    % basic parameter %
    particle_num = 50;    % 粒子群規模
    particle_dim = 2;    % 目標函数的自變量個數
    max_ite = 100;    % 最大迭代次數

    % parameter for disturbance strategy %
    gbest_fitness = 0;    % 使第一次迭代不會出問題 %
    disturb_count = 0;    % 觸發 disturbance strategy 的累積量 %
    disturb_threshold = 5;    % 觸發 disturbance strategy 的 threshold %

    % parameter for stopping criteria %
    gbest_fitness_threshold = 10^(-4);    % 停止迭代條件 %
    stopping_count = 0;    % 觸發 stopping 的累積量 %
    stopping_threshold = 10;    % 觸發 stopping 的 threshold %

    % parameter for inertia weight %
    x_inertia = 0.4;
    A_inertia = 4;
    max_inertia = 0.9;
    min_inertia = 0.4;

    % parameter for updating velocity %
    c1 = 2;
    c2 = 2;
    vmax = 5;    % 使 velocity 的絕對值小於 vmax %

    % parameter for adaptive position update strategy %
    b = 0.3;

    % randomly initialize particle locations and velocities %
    v = 2*rand(particle_num, particle_dim);    % 粒子飛行速度
    x = [0.2+0.6*rand(particle_num, 1) 10^(-1)*rand(particle_num, 1)];    % 粒子所在位置
    pbest_x = x;
    pbest_fitness = 10^(5)*ones(particle_num, 1);    % 先設一個很大的值則第一次迭代時會直接更新

    % 迭代過程 %
    ite = 1;
    while (ite <= max_ite)
        % calculate fitness of particle and update pbest %
        error = zeros(particle_num, 1);
        for i = 1:particle_num
            error(i, :) = optimization_obj_func_Kalman(L, x(i, :), start_ini_frame, ini_frame, Y_DAS, n, Y_delay, look_mic);
            if error(i, :) < pbest_fitness(i, :)
			    pbest_fitness(i, :) = error(i, :);
			    pbest_x(i, :) = x(i, :);
		    end

        end

        % update gbest with disturbance strategy and check gbest 有無變化 %
        gbest_fitness_old = gbest_fitness;
        [gbest_fitness, ii] = min(pbest_fitness);
	    gbest_x = pbest_x(ii, :);

        if gbest_fitness == gbest_fitness_old
            disturb_count = disturb_count + 1;
            stopping_count = stopping_count + 1;
        else
            disturb_count = 0;
            stopping_count = 0;
        end

        if disturb_count == disturb_threshold
            rand_for_disturb = rand;
            Nbest_x = rand_for_disturb*gbest_x + (1-rand_for_disturb)*(gbest_x - pbest_x(randi(particle_num), :));
            Nbest_fitness = optimization_obj_func_Kalman(L, Nbest_x, start_ini_frame, ini_frame, Y_DAS, n, Y_delay, look_mic);
            if Nbest_fitness < gbest_fitness
                gbest_fitness = Nbest_fitness;
                gbest_x = Nbest_x;
                disturb_count = 0;
                stopping_count = 0;
            end

        end

        % check stopping criteria %
        if gbest_fitness < gbest_fitness_threshold || stopping_count == stopping_threshold
            break
        end
        
        % update inertia %
        x_inertia = A_inertia*x_inertia*(1-x_inertia);
        inertia = (max_inertia-min_inertia)*(max_ite-ite)/max_ite + min_inertia*x_inertia;

        % elite and dimensional learning strategies %
        rand_select = randi(particle_num, 4, 1);
        [Cpbest_fitness, index_min] = min(pbest_fitness(rand_select, :));
        Cpbest_x = pbest_x(rand_select(index_min, :), :);

        Fpbest_x = zeros(particle_num, particle_dim);
        for i = 1:particle_num
            if Cpbest_fitness < pbest_fitness(i, :)
                Fpbest_x(i, :) = Cpbest_x;
            else
                Fpbest_x(i, :) = pbest_x(i,:);
            end

        end

        Mpbest_x = mean(pbest_x, 1);

        % update velocity %
        for i = 1:particle_num
		    v(i, :) = inertia*v(i, :) + c1*rand*(Fpbest_x(i, :)-x(i, :)) + c2*rand*(Mpbest_x-x(i,:));

		    for j= 1:particle_dim
			    if v(i, j) > vmax
				    v(i, j) = vmax;
			    elseif v(i, j) < -vmax
				    v(i, j) = -vmax;
                end

            end

        end

        % adaptive position update strategy %
        for i = 1:particle_num
            lambda = exp(error(i, :))/exp(mean(error, 1));
            l = rand*2 - 1;
            if lambda < rand
                x(i, :) = norm(gbest_x-x(i, :))*exp(b*l)*cos(2*pi*l) + gbest_x;
            else
                x(i, :) = x(i, :) + v(i, :);
            end
        end

        % competitive substitution mechanism %
        [W_fitness, index_max] = max(error);
        rand_select = randi(particle_num, 2, 1);
        N_x = gbest_x + rand*(pbest_x(rand_select(1, :), :)-pbest_x(rand_select(2, :), :));
        if optimization_obj_func_Kalman(L, N_x, start_ini_frame, ini_frame, Y_DAS, n, Y_delay, look_mic) < W_fitness
            x(i, :) = N_x;
        end

        %　增加迭代次數 %
        ite = ite + 1;
    end    % end for while iteration

    % 計算最終 weight %
    weight = zeros(L, 1);
    K = zeros(L, 1);    % Kalman gain %
    P = gbest_x(:, 1)*eye(L);    % error covariance matrix %
    R = gbest_x(:, 2);    % measurement noise covariance matrix %
    for FrameNo = start_ini_frame:ini_frame
        % no time update only have measurement update %
        K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
        weight = weight + K*(conj(Y_delay(n, FrameNo, look_mic)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
        P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;
    end

    A(look_mic, :, n) = weight';
end     % end for parfor frequency

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

%  算 NRMSPM %
h_NRMSPM = reshape(h.', [MicNum*points_rir 1]);
aa_NRMSPM = reshape(A_tdomain.', [MicNum*points_rir 1]);
NRMSPM = 20*log10(norm(h_NRMSPM-h_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(h_NRMSPM));

NRMSPM_in = zeros(MicNum, 1);
for i = 1:MicNum
    NRMSPM_in(i, :) = 20*log10(norm(h(i, :).'-h(i, :)*A_tdomain(i, :).'/(A_tdomain(i, :)*A_tdomain(i, :).')*A_tdomain(i, :).')/norm(h(i, :).'));
end

% plot and save estimated RIR and ATF fig %
algorithm = 'ASPSO';
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
ATF = fft(h, points_rir, 2);
ATF_estimated = fft(A_tdomain, points_rir, 2);

figure(2)
plot(h(look_mic, :), 'r');
hold on
plot(A_tdomain(look_mic, :), 'b');
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
xlim([1 points_rir])
legend('ground-truth RIR', 'estimated RIR')
xlabel('time samples')
ylabel('RIR')

fig_filename_str = ['fig\', algorithm, '_', string(reverberation_time), '_RIR.fig'];
fig_filename = join(fig_filename_str, '');
savefig(fig_filename)

figure(3)
subplot(2, 1, 1);
semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF(look_mic, 1:points_rir/2+1))), 'r');
hold on
semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF_estimated(look_mic, 1:points_rir/2+1))), 'b');
hold off
xlim([200 8000])
legend('ground-truth ATF', 'estimated ATF')
xlabel('frequency (Hz)')
ylabel('dB')

subplot(2, 1, 2);
semilogx(linspace(0, fs/2, points_rir/2+1), unwrap(angle(ATF(look_mic, 1:points_rir/2+1))), 'r');
hold on
semilogx(linspace(0, fs/2, points_rir/2+1), unwrap(angle(ATF_estimated(look_mic, 1:points_rir/2+1))), 'b');
hold off
xlim([200 8000])
legend('ground-truth ATF', 'estimated ATF')
xlabel('frequency (Hz)')
ylabel('phase (radius)')

fig_filename_str = ['fig\', algorithm, '_', string(reverberation_time), '_ATF.fig'];
fig_filename = join(fig_filename_str, '');
savefig(fig_filename)

% save A_tdomain %
A_filename_str = ['A_tdomain\', algorithm, '_', string(reverberation_time), '_A_tdomain.mat'];
A_filename = join(A_filename_str, '');
save(A_filename, 'A_tdomain')

fprintf('done\n')
