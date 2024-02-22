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
points_rir = 12288;                                      % Number of rir points (需比 reverberation time 還長)
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 1;                                           % Disable high-pass filter

%% generate ground-truth RIR (h) %%
% 產生 RIR 和存.mat 檔 %
% h = rir_generator(c, fs, MicPos, SorPos, room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
% rir_filename_str = ['h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
% rir_filemane = join(rir_filename_str, '');
% save(rir_filemane, 'h')

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

%% WPE %%
% % do wpe %
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
% y_wpe = y_wpe.';
% 
% % 存 wpe mat %
% y_wpe_filename_str = ['y_wpe-', string(reverberation_time), '.mat'];
% y_wpe_filename = join(y_wpe_filename_str, '');
% save(y_wpe_filename, 'y_wpe')
% 
% % 存 wpe wav %
% ratio_y_wpe = 0.8 / max(abs(y_wpe(look_mic, :))) ;
% y_filemane_str = ['wav\y_wpe_all-', string(reverberation_time), '.wav'];
% y_filemane = join(y_filemane_str, '');
% audiowrite(y_filemane, y_wpe(look_mic, :)*ratio_y_wpe, fs)

% load y_wpe %
y_wpe_filename_str = ['y_wpe-', string(reverberation_time), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

% ratio_y = 0.8 / max(abs(y_nodelay(look_mic, :))) ;
% ratio_y_wpe = 0.8 / max(abs(y_wpe(look_mic, :))) ;
% % 畫 y_wpe time plot %
% figure(5);
% plot(y_nodelay(look_mic, :), 'r');
% hold on
% plot(y_wpe(look_mic, :), 'b');
% hold off
% ylim([-1.1 1.1])
% title('y\_wpe')
% xlabel('points')
% ylabel('magnitude')
% legend('y\_nodelay', 'y\_wpe')
% shg

%% DAS beamformer %%
% 算 mic 與 source 之距離 %
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((SorPos - MicPos(i, :)).^2));
end

% 算 a %
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

%% PSO %%
start_ini_frame = ceil((12*fs - NFFT)/hopsize) + 1;    % wpe 第12秒開始穩定
ini_frame = NumOfFrame;
A = zeros(MicNum, L, frequency);

tic
parfor n =1:frequency
    % basic parameter %
    particle_num = 50;    % 粒子群規模
    particle_dim = 2;    % 目標函数的自變量個數
    max_ite = 100;    % 最大迭代次數
    
    % parameter for stopping criteria %
    gbest_fitness_threshold = 10^(-4);    % 停止迭代條件 %
    gbest_fitness = 0;    % 使第一次迭代不會出問題 %
    stopping_count = 0;    % 觸發 stopping 的累積量 %
    stopping_threshold = 10;    % 觸發 stopping 的 threshold %

    % parameter for updating velocity %
    inertia_factor = 0.6;    % 惯性因子
    c1 = 2;    % 每個粒子的個體學習因子，加速度常數
    c2 = 2;    % 每個粒子的社會學習因子，加速度常數
    vmax = 5;    % 粒子的最大飛行速度

    % randomly initialize particle locations and velocities %
    v = 2*rand(particle_num, particle_dim);    % 粒子飛行速度
    x = [0.2+0.6*rand(particle_num, 1) 10^(-1)*rand(particle_num, 1)];    % 粒子所在位置
    pbest_x = x;
    pbest_fitness = 10^(5)*ones(particle_num, 1);    % 先設一個很大的值則第一次迭代時會直接更新

    ite=1;
    while (ite <= max_ite)
        % calculate fitness of particle and update pbest %
        error = zeros(particle_num, 1);
        for i = 1:particle_num
            error(i, :) = obj_func_Kalman(L, x(i, :), start_ini_frame, ini_frame, Y_DAS, n, Y_delay, look_mic);
            if error(i, :) < pbest_fitness(i, :)
			    pbest_fitness(i, :) = error(i, :);
			    pbest_x(i,:) = x(i,:);
		    end

        end

        % update gbest %
        gbest_fitness_old = gbest_fitness;
        [gbest_fitness, ii] = min(pbest_fitness);
	    gbest_x = pbest_x(ii, :);

         % check gbest 有無變化 %
        if gbest_fitness == gbest_fitness_old
            stopping_count = stopping_count + 1;
        else
            stopping_count = 0;
        end

        % check stopping criteria %
        if gbest_fitness < gbest_fitness_threshold || stopping_count == stopping_threshold
            break
        end

        % update velocity and location %
        for i = 1:particle_num
		    v(i, :) = inertia_factor*v(i, :) + c1*rand*(pbest_x(i, :)-x(i, :)) + c2*rand*(gbest_x-x(i,:));
		    for j= 1:particle_dim
			    if v(i, j) > vmax
				    v(i, j) = vmax;
			    elseif v(i, j) < -vmax
				    v(i, j) = -vmax;
                end

            end

		    x(i,:) = x(i,:) + v(i,:);
        end
       
        %　增加迭代次數 %
        ite = ite + 1;
    end    % end for while iteration

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

A_dB = mag2db(abs(A));
% 畫 A frequency plot %
figure(6);
mesh(L_vector, freqs_vector, squeeze(A_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% A 轉回時域 (A_tdomain) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
ratio_A_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_tdomain(i, :) = max(abs(h(i, :)))/max(abs(A_tdomain(i, :)));
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
ATF = fft(h, points_rir, 2);
ATF_estimated = fft(A_tdomain(:, hopsize*(osfac-1)+1:end), points_rir, 2);
sum_norm = 0;
for i  = look_mic:look_mic
    norm_ATF = norm(ATF(i, :) - ATF_estimated(i, :));
    sum_norm = sum_norm + norm_ATF;
end

ME = sum_norm;

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
