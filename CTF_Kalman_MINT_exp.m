clc; clear;
close all;

%% RIR parameter %%
SorNum = 1;                                              % source number
MicNum = 6;                                             % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 48000;                                              % Sample frequency (samples/s)
Ts = 1/fs;                                               % Sample period (s)
points_rir = 8192;                                      % Number of rir points (需比 reverberation time 還長)

%% window parameter %%
NFFT = 1024;
hopsize = 256;

win = hamming(NFFT);
osfac = round(NFFT/hopsize);

frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);
L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);

%% 讀音檔 or 產生 white noise source (source) %%
Second = 8;
SorLen =  Second*fs;

% load source %
[source_transpose, fs] = audioread('exp_wav\loudspeasker_1.wav', [1, SorLen]);    % speech source
source = source_transpose.';

%% compute source signal for frequency (S) %%
% source 轉頻域 %
source_transpose = source.';
[S, ~, S_t_vector] = stft(source_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');

NumOfFrame = size(S_t_vector, 1);
NumOfFrame_vector = 1:1:NumOfFrame;
S_dB = mag2db(abs(S));
% % 畫 source frequency plot %
% figure(1);
% mesh(NumOfFrame_vector, freqs_vector, S_dB)
% colorbar
% view(2)
% xlim([1 NumOfFrame])
% title('S')
% xlabel('frame')
% ylabel('frequency(Hz)')
% shg

%% RIR mix source 先在時域上 convolution 再做 stft (y and Y) %%
for i = 1:MicNum
    y_nodelay_str = ['exp_wav\', string(i), 'th.wav'];
    y_nodelay_filename = join(y_nodelay_str, '');
    [y_nodelay(i, :), fs] = audioread( y_nodelay_filename, [1, SorLen]); 
end
extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y_delay = zeros(MicNum, SorLen);
y_delay(:, extra_delay_y+1:end) = y_nodelay(:, 1:SorLen-extra_delay_y);


% y 轉頻域 %
y_delay_transpose = y_delay.';
[Y_delay, ~, ~] = stft(y_delay_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');


%% predict CTF with Kalman filter (A) %%
start_ini_frame = L;
ini_frame = NumOfFrame;

% Kalman filter %
A = zeros(MicNum, L, frequency);
tic
parfor i = 1:MicNum
    for n =1:frequency
        weight = zeros(L, 1);
        P = 0.5*eye(L);    % error covariance matrix %
        K = zeros(L, 1);    % Kalman gain %
        R = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_ini_frame:ini_frame
            % no time update only have measurement update %
            K = P*flip(S(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(S(n, FrameNo-L+1:FrameNo)))*P*flip(S(n, FrameNo-L+1:FrameNo).') + R);
            weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(S(n, FrameNo-L+1:FrameNo)))*weight);
            P = P - K*conj(flip(S(n, FrameNo-L+1:FrameNo)))*P;
        end
    
        A(i, :, n) = weight';
    end

end
toc

look_mic = 3;
A_dB = mag2db(abs(A));
% 畫 A frequency plot %
figure(2);
mesh(L_vector, freqs_vector, squeeze(A_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% TF estimate %%
TF =  tfestimate(source, y_nodelay.', hamming(points_rir), points_rir*0.75, points_rir);
TF = TF.';
TF = [TF flip(conj(TF(:, 2:end-1)))];
tf = ifft(TF, points_rir, 2,'symmetric');

e_total = sum(tf.^2, 2);
t60 = zeros(MicNum, 1);
for i = 1:MicNum
    edc = zeros(points_rir, 1);
    for n = 1:points_rir
        edc(n, :) = 10*log10(sum(tf(i, n:end).^2)/e_total(i, :));
    end
    
    edc = edc + 60;
    [~, t60_point] = min(abs(edc));
    t60(i, :) = t60_point/fs;
end

figure(3)
plot(tf(look_mic, :));
title('tfestimate')
xlabel('points')
ylabel('amplitude')
shg

%% A 轉回時域 (A_tdomain) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
ratio_A_tdomain = zeros(MicNum, 1);
for i = 1:MicNum
    ratio_A_tdomain(i, :) = max(abs(tf(i, :)))/max(abs(A_tdomain(i, hopsize*(osfac-1)+1:end)));
end

A_tdomain = A_tdomain.*ratio_A_tdomain;

A_tdomain_forplot = A_tdomain(look_mic, hopsize*(osfac-1)+1:end);
% 畫 A_tdomain time plot
figure(4)
plot(tf(look_mic, :), 'r');
hold on
plot(A_tdomain_forplot, 'b');
hold off
legend('MTF', 'CTF')
xlabel('points')
ylabel('amplitude')
shg
