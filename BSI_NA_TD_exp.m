clc; clear;
close all;

tic

%% RIR parameter %%
SorNum = 1;                                              % source number
MicNum = 6;                                              % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 48000;                                              % Sample frequency (samples/s)
Ts = 1/fs;                                               % Sample period (s)
points_rir = 8192;                                       % Number of rir points (需比 reverberation time 還長)
look_mic = 6;

%% read source 音檔 (source) %%
Second = 8;
SorLen =  Second*fs;

% load source %
[source_transpose, Fs] = audioread('wav_exp\loudspeaker_1.wav', [1, SorLen]);    % speech source
source = source_transpose.';

%% read mic 音檔 (x) %%
% load x %
x = zeros(MicNum, SorLen);
for i = 1:MicNum
    x_str = ['wav_exp\', string(i-1), 'th.wav'];
    x_filename = join(x_str, '');
%     [x(i, :), fs] = audioread( x_filename, [1, SorLen]);    % 第一顆喇叭
    [x(i, :), fs] = audioread( x_filename, [30*fs+1, 30*fs+SorLen]);    % 第三顆喇叭
end

%% TF estimate (tf) %%
TF =  tfestimate(source, x.', hamming(points_rir), points_rir*0.75, points_rir);
TF = TF.';
TF = [TF flip(conj(TF(:, 2:end-1)))];
tf = ifft(TF, points_rir, 2,'symmetric');
h = tf;

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
for n = L+1:SorLen
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
    ratio_aa(i, :) = max(abs(h(i, :)))/max(abs(aa(i, :)));
end

aa = aa.*ratio_aa;

look_mic = 1;
figure(3)
plot(h(look_mic, :), 'r');
hold on
plot(-aa(look_mic, :), 'b');
hold off
xlim([1 points_rir])
legend('ground-truth RIR', 'estimated RIR')
xlabel('time samples')
ylabel('amplitude')
shg

h_NRMSPM = reshape(h.', [MicNum*L 1]);
aa_NRMSPM = reshape(aa.', [MicNum*L 1]);
NRMSPM = 20*log(norm(h_NRMSPM-h_NRMSPM.'*aa_NRMSPM/(aa_NRMSPM.'*aa_NRMSPM)*aa_NRMSPM)/norm(h_NRMSPM));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reshape and rescale h_hat %
h_hat = reshape(h_hat, [size(h_hat, 1)/MicNum MicNum]).';
ratio_h_hat = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_h_hat(i, :) = max(abs(h(i, :)))/max(abs(h_hat(i, :)));
end

h_hat = h_hat.*ratio_h_hat;

% 畫圖看結果 %
figure(4)
plot(cost_fun(L+1:end, :).');
xlabel('update times')
title('cost function')

look_mic = 1;
figure(5)
plot(h(look_mic, :), 'r');
hold on
plot(-h_hat(look_mic, :), 'b');
hold off
xlim([1 points_rir])
legend('ground-truth RIR', 'estimated RIR')
xlabel('time samples')
ylabel('amplitude')
shg

% ME %
ATF = fft(h, points_rir, 2);
ATF_estimated = fft(h_hat, points_rir, 2);
sum_norm = 0;
for i  = 1:MicNum
    norm_ATF = norm(ATF(i, :) - ATF_estimated(i, :));
    sum_norm = sum_norm + norm_ATF;
end

ME = sum_norm/MicNum;

% NRMSPM %
h_NRMSPM = reshape(h.', [MicNum*L 1]);
h_hat_NRMSPM = reshape(h_hat.', [MicNum*L 1]);
NRMSPM = 20*log(norm(h_NRMSPM-h_NRMSPM.'*h_hat_NRMSPM/(h_hat_NRMSPM.'*h_hat_NRMSPM)*h_hat_NRMSPM)/norm(h_NRMSPM));

toc
