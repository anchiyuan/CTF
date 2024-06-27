clc; clear;
close all;

%% RIR parameter %%
SorNum = 1;
c = 343;
Fs = 48000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 16000;           % 欲 resample 成的取樣頻率
MicNum = 3;           % 跑不了太多通道因為記憶體不夠
points_rir = 2048;    % 自行設定想要輸出的 RIR 長度
date = '0607';
look_mic = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read source 音檔 (source) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Second = 28;             % 使用時間長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*Fs;
sorLen =  Second*fs;

% load source %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source_str = ['wav_exp\', date,'\white_noise_30s.wav'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source_filename = join(source_str, '');
[source_transpose, ~] = audioread(source_filename, [1, SorLen]);
source = source_transpose.';

% resample %
source = resample(source, 1, Fs/fs);

% load source %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
speaker_str = ['wav_exp\', date,'\14.wav'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
speaker_filename = join(speaker_str, '');
[speaker_transpose, ~] = audioread(speaker_filename, [1, SorLen]);
speaker = speaker_transpose.';

% resample %
speaker = resample(speaker, 1, Fs/fs);

%% read mic 音檔 (x) %%
% load x %
x = zeros(MicNum, SorLen);
for i = 1:MicNum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_str = ['wav_exp\', date,'\', string(i+6), '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_filename = join(x_str, '');
    [x(i, :), ~] = audioread( x_filename, [1, SorLen]);
end
% resample %
x = resample(x, 1, Fs/fs, Dimension=2);

%% TF estimate (tf) %%
TF =  tfestimate(speaker, x.', hamming(points_rir), points_rir*0.75, points_rir);
TF = TF.';
TF = [TF flip(conj(TF(:, 2:end-1)))];
tf = ifft(TF, points_rir, 2,'symmetric');

figure(1)
plot(tf(look_mic, :));
title('tfestimate')
xlabel('points')
ylabel('amplitude')
shg

%% MCN %%
% basic parameters %
L = points_rir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 0.999;
step_size = 0.999;
diaload_delta_h = 0;
stop_criteria = 1e-8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize R %
R_ini_scalar = 0;
for i = 1:MicNum
    R_ini_scalar = R_ini_scalar + flip(x(i, 1:L))*flip(x(i, 1:L)).'/L;
end

R = R_ini_scalar*eye(MicNum*L);

% initialize h_hat %
h_hat = zeros(MicNum*L, 1);
for i = 1:MicNum
    h_hat((i-1)*L+1) = 1;
end

h_hat = h_hat/norm(h_hat);

% initialize cost function %
cost_fun = zeros(SorLen, 1);

% iteration process %
for n = 7858:SorLen
    % construct R %
    R_temp = zeros(MicNum*L, MicNum*L);
    dia_sum = zeros(L, L);
    for i = 1:MicNum
        for j = 1:MicNum
            R_temp((i-1)*L+1:i*L, (j-1)*L+1:j*L) = -(flip(x(j, n-L+1:n)).'*flip(x(i, n-L+1:n)));
            if i == j
                dia_sum = dia_sum + flip(x(j, n-L+1:n)).'*flip(x(i, n-L+1:n));
            end

        end

    end

    for i = 1:MicNum
        R_temp((i-1)*L+1:i*L, (i-1)*L+1:i*L) = R_temp((i-1)*L+1:i*L, (i-1)*L+1:i*L) + dia_sum;
    end

    R = lambda*R + R_temp;

    % construct cost function %
    for i = 1:MicNum-1
        for j = i+1:MicNum
            cost_fun(n, :) = cost_fun(n, :) + (flip(x(i, n-L+1:n))*h_hat((j-1)*L+1:j*L, :) - flip(x(j, n-L+1:n))*h_hat((i-1)*L+1:i*L, :))^2;
        end

    end

    % compute new h_hat %
    delta_h = inv(R-2*h_hat*(h_hat.')*R-2*R*h_hat*(h_hat.')+diaload_delta_h*eye(MicNum*L))*(R_temp*h_hat-cost_fun(n, :)*h_hat);    % diaload_delta_h*eye(MicNum*L)
    h_hat = (h_hat-step_size*delta_h)/norm(h_hat-step_size*delta_h);

    % stop criteria %
    cost_mean = mean(cost_fun(n-9:n, :));
    if cost_mean < stop_criteria
        break;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 中途看波型跟NRMSPM %
aa = reshape(h_hat, [size(h_hat, 1)/MicNum MicNum]).';
ratio_aa = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_aa(i, :) = max(abs(tf(i, :)))/max(abs(aa(i, :)));
end

aa = aa.*ratio_aa;

figure(2)
plot(tf(look_mic, :), 'r');
hold on
plot(-aa(look_mic, :), 'b');
hold off
xlim([1 points_rir])
legend('tfestimate', 'BSI')
title('RIR')
xlabel('time samples')
ylabel('amplitude')
shg

tf_NRMSPM = reshape(tf.', [MicNum*L 1]);
aa_NRMSPM = reshape(aa.', [MicNum*L 1]);
NRMSPM_mid = 20*log10(norm(tf_NRMSPM-tf_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(tf_NRMSPM));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 畫圖看結果 %%
% reshape and rescale h_hat %
h_hat = reshape(h_hat, [size(h_hat, 1)/MicNum MicNum]).';
ratio_h_hat = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_h_hat(i, :) = max(abs(tf(i, :)))/max(abs(h_hat(i, :)));
end

h_hat = h_hat.*ratio_h_hat;

A_filename = 'A_tdomain_exp\BSI_NA_TD_h_hat.m';
save(A_filename, 'h_hat')

% cost function 圖 %
figure(3)
plot(cost_fun(L+1:n, :).');
xlabel('update times')
title('cost function')

fig_filename = 'fig_exp\BSI_NA_TD_costfun.fig';
savefig(fig_filename)

% RIR 比較圖 %
figure(4)
plot(tf(look_mic, :), 'r');
hold on
plot(-h_hat(look_mic, :), 'b');
hold off
xlim([1 points_rir])
legend('tfestimate', 'BSI')
title('RIR')
xlabel('time samples')
ylabel('amplitude')
shg

fig_filename = 'fig_exp\BSI_NA_TD_RIR.fig';
savefig(fig_filename)

%% NRMSPM %%
tf_NRMSPM = reshape(tf.', [MicNum*points_rir 1]);
h_hat_NRMSPM = reshape(h_hat.', [MicNum*points_rir 1]);
NRMSPM = 20*log10(norm(tf_NRMSPM-tf_NRMSPM.'*h_hat_NRMSPM/(h_hat_NRMSPM.'*h_hat_NRMSPM)*h_hat_NRMSPM)/norm(tf_NRMSPM));

NRMSPM_in = zeros(MicNum, 1);
for i = 1:MicNum
    NRMSPM_in(i, :) = 20*log10(norm(tf(i, :).'-tf(i, :)*h_hat(i, :).'/(h_hat(i, :)*h_hat(i, :).')*h_hat(i, :).')/norm(tf(i, :).'));
end

%% 檢查 direct sound 有無 match %%
[~, argmax_tf] = max(abs(tf.'));
[~, argmax_h_hat] = max(abs(h_hat.'));

%% ATF magnitude and phase plot %%
ATF = fft(tf, points_rir, 2);
ATF_estimated = fft(h_hat, points_rir, 2);

figure(5)
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

fig_filename = 'fig_exp\BSI_NA_TD_ATF.fig';
savefig(fig_filename)

fprintf('done\n')
