clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')

%% RIR parameter %%
SorNum = 1;
c = 343;
Fs = 48000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 16000;           % 欲 resample 成的取樣頻率
MicNum = 8;           % 實驗麥克風數量
points_rir = 2048;    % 自行設定想要輸出的 RIR 長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% window parameter %%
NFFT = 1024;
hopsize = 256;

win = hamming(NFFT);
osfac = round(NFFT/hopsize);

frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);
freqs_vector = linspace(0, fs/2, frequency);

%% read source 音檔 (source) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Second = 28;             % 使用時間長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*Fs;
sorLen =  Second*fs;

% load source %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[source_transpose, ~] = audioread('wav_exp\white_noise_30s.wav', [1, SorLen]);    % speech source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source = source_transpose.';

% resample %
source = resample(source, 1, Fs/fs);

%% source 做 stft (S) %%
[S, ~, ~] = stft(source.', fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

NumOfFrame = size(S, 2);

%% read mic 音檔再做 stft (y_nodelay, y_delay and Y_delay) %%
% load y_nodelay %
y_nodelay = zeros(MicNum, SorLen);
for i = 1:MicNum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_nodelay_str = ['wav_exp\', string(i), '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_nodelay_filename = join(y_nodelay_str, '');
    [y_nodelay(i, :), ~] = audioread( y_nodelay_filename, [1, SorLen]);
end

% resample %
y_nodelay = resample(y_nodelay, 1, Fs/fs, Dimension=2);

% delay y_nodelay to get y_delay %
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, sorLen);
y_delay(:, extra_delay_y+1:end) = y_nodelay(:, 1:sorLen-extra_delay_y);

% y_delay 轉頻域 to get Y_delay %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPE (y_wpe) %%
% % do wpe %
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
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

%% TDOA localization %%
% GCC-PHAT for delay estimation %
delay = zeros(MicNum-1, 1);
difference = zeros(MicNum-1, 1);
for i = 1:MicNum-1
    delay(i, :) = gccphat(y_nodelay(i+1,:).', y_nodelay(1,:).', fs);
    difference(i, :) = delay(i, :)*c;
end

% mics position with repect to reference point %
mic_x = [ 0 ; 91.8 ; 91.8 ;    0 ;    0 ; 91.8 ; 91.8 ;    0 ]./100;
mic_y = [ 0 ;    0 ; 90.4 ; 90.6 ;    0 ;    0 ; 90.4 ; 90.6 ]./100;
mic_z = [ 0 ;    0 ;    0 ;    0 ; 80.4 ; 79.8 ;   80 ;   80 ]./100;

micpos = [mic_x, mic_y, mic_z,];

% generate parameters matrix %
A = 2*[micpos(2:end, :), difference];
b = micpos(2:end, 1).^2 + micpos(2:end, 2).^2 + micpos(2:end, 3).^2 - difference.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = 4*sigma^2*(ones(MicNum-1, 1)*ones(MicNum-1, 1).' + eye(MicNum-1))*difference*difference.';
dia_load_psi = 1;
invpsi = inv(psi+dia_load_psi*eye(MicNum-1)); 
P = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 -1];

% find possible root %
syms Z
theta = pinv(A.'*invpsi*A+Z*P )*(A.'*invpsi*b);
I = theta.'*P*theta;
[num, den] = numden(I);
poly = sym2poly(num);
roots = roots(poly);
roots_real = roots(imag(roots)==0, :);    % 取只有實數的根
if isempty(roots_real)
    roots_real = [roots_real; 0];
end

% find actual root from posible root %
if size(roots_real, 1) > 1
    theta = zeros(4, size(roots_real, 1));
    costfun = zeros(size(roots_real, 1), 1);
    for i = 1:size(roots_real, 1)
        theta(:, i) = pinv(A.'*invpsi*A+roots_real(i, :)*P )*(A.'*invpsi*b);
        costfun(i, :) = (A*theta(:, i)-b).'*invpsi*(A*theta(:, i)-b);
    end

    [~, min_index] = min(costfun);
    sorpos_estimation = theta(1:3, min_index).';

else
    theta = pinv(A.'*invpsi*A+roots_real*P )*(A.'*invpsi*b);
    sorpos_estimation = theta(1:3, :).';

end

%% TF estimate (tf and t60) %%
TF =  tfestimate(source, y_nodelay.', hamming(points_rir), points_rir*0.75, points_rir);
TF = TF.';
TF = [TF flip(conj(TF(:, 2:end-1)))];
tf = ifft(TF, points_rir, 2,'symmetric');

e_total = sum(tf.^2, 2);
t60 = zeros(MicNum, 1);
for i = 1:MicNum
    edc = zeros(points_rir, 1);
    for j = 1:points_rir
        edc(j, :) = 10*log10(sum(tf(i, j:end).^2)/e_total(i, :));
    end
    
    edc = edc + 60;
    [~, t60_point] = min(abs(edc));
    t60(i, :) = t60_point/fs;
end

look_mic = 1;
figure(1)
plot(tf(look_mic, :));
title('tfestimate')
xlabel('points')
ylabel('amplitude')
shg

%% ASPSO %%
y_wpe_transpose = y_wpe.';
[Y_wpe, ~, ~] = stft(y_wpe_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

% basic parameter %
particle_num = 50;    % 粒子群規模
particle_dim = 3;    % 目標函数的自變量個數
max_ite = 100;    % 最大迭代次數

% parameter for disturbance strategy %
gbest_fitness = 0;    % 使第一次迭代不會出問題 %
disturb_count = 0;    % 觸發 disturbance strategy 的累積量 %
disturb_threshold = 5;    % 觸發 disturbance strategy 的 threshold %

% parameter for stopping criteria %
gbest_fitness_threshold = -10;    % 停止迭代條件 %
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
x = [sorpos_estimation(:, 1)+0.2*(rand(particle_num, 1)*2-1), sorpos_estimation(:, 2)+0.2*(rand(particle_num, 1)*2-1), sorpos_estimation(:, 3)+0.2*(rand(particle_num, 1)*2-1)];    % 粒子所在位置
pbest_x = x;
pbest_fitness = 10^(5)*ones(particle_num, 1);    % 先設一個很大的值則第一次迭代時會直接更新

% 迭代過程 %
ite = 1;
while (ite <= max_ite)
    % calculate fitness of particle and update pbest %
    error = zeros(particle_num, 1);
    for i = 1:particle_num
        error(i, :) = obj_func_exp(MicNum, SorNum, x(i, :), micpos, frequency, freqs_vector, c, NumOfFrame, Y_wpe, L, Y_delay, points_rir, NFFT, hopsize, win, fs, osfac, tf);
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
        Nbest_fitness = obj_func_exp(MicNum, SorNum, Nbest_x, micpos, frequency, freqs_vector, c, NumOfFrame, Y_wpe, L, Y_delay, points_rir, NFFT, hopsize, win, fs, osfac, tf);
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
    if obj_func_exp(MicNum, SorNum, N_x, micpos, frequency, freqs_vector, c, NumOfFrame, Y_wpe, L, Y_delay, points_rir, NFFT, hopsize, win, fs, osfac, tf) < W_fitness
        x(i, :) = N_x;
    end

    %　增加迭代次數 %
    ite = ite + 1;
end    % end for while iteration

%% A 轉回時域畫圖 (A_tdomain) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
A_tdomain = A_tdomain(:, hopsize*(osfac-1)+1:end);
ratio_A_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_tdomain(i, :) = max(abs(tf(i, :)))/max(abs(A_tdomain(i, :)));
end

A_tdomain = A_tdomain.*ratio_A_tdomain;

look_mic = 1;
% 畫 A_tdomain time plot
figure(2)
plot(tf(look_mic, :), 'r');
hold on
plot(-A_tdomain(look_mic, :), 'b');
hold off
legend('tfestimate', 'CTF')
title('RIR')
xlabel('points')
ylabel('amplitude')
shg

%% 算 NRMSPM %%
tf_NRMSPM = reshape(tf.', [MicNum*points_rir 1]);
A_tdomain_NRMSPM = reshape(A_tdomain.', [MicNum*points_rir 1]);
NRMSPM = 20*log10(norm(tf_NRMSPM-tf_NRMSPM.'*A_tdomain_NRMSPM/(A_tdomain_NRMSPM.'*A_tdomain_NRMSPM)*A_tdomain_NRMSPM)/norm(tf_NRMSPM));

NRMSPM_in = zeros(MicNum, 1);
for i = 1:MicNum
    NRMSPM_in(i, :) = 20*log10(norm(tf(i, :).'-tf(i, :)*A_tdomain(i, :).'/(A_tdomain(i, :)*A_tdomain(i, :).')*A_tdomain(i, :).')/norm(tf(i, :).'));
end

%% 檢查 direct sound 有無 match %%
[~, argmax_tf] = max(abs(tf.'));
[~, argmax_A_tdomain] = max(abs(A_tdomain.'));

%% 檢查 ATF 有無 match %%
ATF = fft(tf, points_rir, 2);
ATF_estimated = fft(A_tdomain, points_rir, 2);

figure(3)
semilogx(linspace(0, fs/2, points_rir/2+1), abs(ATF(look_mic, 1:points_rir/2+1)), 'r');
hold on
semilogx(linspace(0, fs/2, points_rir/2+1), abs(ATF_estimated(look_mic, 1:points_rir/2+1)), 'b');
hold off
legend('tfestimate', 'CTF')
title('ATF')
xlabel('frequency')
ylabel('magnitude')
shg

fprintf('done\n')
