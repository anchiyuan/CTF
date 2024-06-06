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
MicNum_TDOA = 8;      % TDOA麥克風數量
MicNum = 6;           % UL麥克風數量
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
% load y_TDOA %
y_TDOA = zeros(MicNum_TDOA, SorLen);
for i = 1:MicNum_TDOA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_TDOA_str = ['wav_exp\', string(i), '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_TDOA_filename = join(y_TDOA_str, '');
    [y_TDOA(i, :), ~] = audioread(y_TDOA_filename, [1, SorLen]);
end

% load y_nodelay %
y_nodelay = zeros(MicNum, SorLen);
for i = 1:MicNum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_nodelay_str = ['wav_exp\', string(i+MicNum_TDOA), '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_nodelay_filename = join(y_nodelay_str, '');
    [y_nodelay(i, :), ~] = audioread( y_nodelay_filename, [1, SorLen]);
end

% resample %
y_TDOA = resample(y_TDOA, 1, Fs/fs, Dimension=2);
y_nodelay = resample(y_nodelay, 1, Fs/fs, Dimension=2);

% delay y_nodelay to get y_delay %
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, sorLen);
y_delay(:, extra_delay_y+1:end) = y_nodelay(:, 1:sorLen-extra_delay_y);

% y_delay 轉頻域 to get Y_delay %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

%% WPE (y_wpe) %%
% do wpe %
y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
y_wpe = y_wpe.';

% 存 wpe mat %
y_wpe_str = ['y_exp\y_wpe_', string(fs),'.mat'];
y_wpe_filename = join(y_wpe_str, '');
save(y_wpe_filename, 'y_wpe')

% load y_wpe %
y_wpe_str = ['y_exp\y_wpe_', string(fs),'.mat'];
y_wpe_filename = join(y_wpe_str, '');
load(y_wpe_filename);

%% TDOA localization %%
% GCC-PHAT for delay estimation %
delay = zeros(MicNum_TDOA-1, 1);
difference = zeros(MicNum_TDOA-1, 1);
for i = 1:MicNum_TDOA-1
    delay(i, :) = gccphat(y_TDOA(i+1,:).', y_TDOA(1,:).', fs);
    difference(i, :) = delay(i, :)*c;
end

% mics position with repect to reference point %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mic_x = [ 0 ;   92 ;   92 ;    0 ;    0 ;   92 ;   92 ;    0 ]./100;
mic_y = [ 0 ;    0 ; 90.8 ;   90 ;    0 ;    0 ; 90.8 ;   90 ]./100;
mic_z = [ 0 ;    0 ;    0 ;    0 ; 80.1 ;   80 ; 80.4 ; 79.9 ]./100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
micpos = [mic_x, mic_y, mic_z,];

% generate parameters matrix %
A = 2*[micpos(2:end, :), difference];
b = micpos(2:end, 1).^2 + micpos(2:end, 2).^2 + micpos(2:end, 3).^2 - difference.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = 4*sigma^2*(ones(MicNum_TDOA-1, 1)*ones(MicNum_TDOA-1, 1).' + eye(MicNum_TDOA-1))*difference*difference.';
dia_load_psi = 1;
invpsi = inv(psi+dia_load_psi*eye(MicNum_TDOA-1)); 
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



%% Ryy %%
% y_nodelay_transpose = y_nodelay.';
% [Y_nodelay, ~, ~] = stft(y_nodelay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
% 
% 
% Ryy = zeros(MicNum, MicNum, frequency);
% for FrameNo= 1:NumOfFrame
%     for n = 1: frequency
%          Ryy(:, :, n) = Ryy(:,:,n) + squeeze(Y_nodelay(n, FrameNo, :))*squeeze(Y_nodelay(n, FrameNo, :))';
%     end  
% 
% end
% 
% Ryy = Ryy/NumOfFrame;

%% DAS or MPDR beamformer (Y_DAS) %%
% 算 mic 與 source 之距離 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MicStart = [0.076, -0.602, 0];
spacing = 0.07;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ULA_pos = zeros(MicNum, 3);
for i = 1:MicNum
    ULA_pos(i, :) = [MicStart(1, 1)+(i-1)*spacing, MicStart(1, 2), MicStart(1, 3)];
end
 
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((sorpos_estimation - ULA_pos(i, :)).^2));
end

% 算 MPDR weight %
% a = zeros(MicNum, SorNum, frequency);
% w = zeros(MicNum, SorNum, frequency);
% dia_load_MPDR = 10^(-4);
% for n = 1:frequency
%     omega = 2*pi*freqs_vector(n);
%     a(:, :, n) = exp(-1j*omega/c*distance)./distance;
%     w(:, :, n) = inv(Ryy(:, :, n)+eye(MicNum)*dia_load_MPDR)*a(:, :, n)/(a(:, :, n)'*inv(Ryy(:, :, n)+eye(MicNum)*dia_load_MPDR)*a(:, :, n));
% 
% end

% DAS weight %
a = zeros(MicNum, SorNum, frequency);
for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, :, n) = exp(-1j*omega/c*distance)./distance;
end

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

%% predict CTF with Kalman filter (A) %%
start_ini_frame = L;
ini_frame = NumOfFrame;
A = zeros(MicNum, L, frequency);

% Kalman stationary filter %
tic
parfor i = 1:MicNum
    for n =1:frequency
        weight = zeros(L, 1);
        P = 0.5*eye(L);    % error covariance matrix %
        K = zeros(L, 1);    % Kalman gain %
        R = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_ini_frame:ini_frame
            % no time update only have measurement update %
            K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
            weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
            P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;
        end
    
        A(i, :, n) = weight';
    end

end
toc

% Kalman nonstationary filter %
% tic
% parfor i = 1:MicNum
%     for n =1:frequency
%         weight = zeros(L, 1);
%         P = 0.5*eye(L);    % error covariance matrix %
%         K = zeros(L, 1);    % Kalman gain %
%         Q = 10^(-4)*eye(L);    % process noise covariance matrix %
%         R = 10^(-3);    % measurement noise covariance matrix %
%         for FrameNo = start_ini_frame:ini_frame
%             % time update %
%             P = P + Q;
% 
%             % measurement update %
%             K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
%             weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
%             P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;
% 
%         end
% 
%         A(i, :, n) = weight';
%     end
% 
% end
% toc

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

%% 畫圖檢查 ATF 有無 match %%
ATF = fft(tf, points_rir, 2);
ATF_estimated = fft(A_tdomain, points_rir, 2);

figure(3)
semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF(look_mic, 1:points_rir/2+1))), 'r');
hold on
semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF_estimated(look_mic, 1:points_rir/2+1))), 'b');
hold off
legend('tfestimate', 'CTF')
title('ATF')
xlabel('frequency')
ylabel('dB')
shg

%% 畫圖檢查 S 和 Y_DAS 有無 match %%
S_dB = mag2db(abs(S));
figure(4);
mesh(1:1:NumOfFrame, freqs_vector, S_dB)
colorbar
view(2)
title('S')
xlabel('frame')
ylabel('frequency(Hz)')
shg

Y_DAS_dB = mag2db(abs(Y_DAS));
figure(5);
mesh(1:1:NumOfFrame, freqs_vector, Y_DAS_dB)
colorbar
view(2)
title('DAS')
xlabel('frame')
ylabel('frequency(Hz)')
shg

fprintf('done\n')
