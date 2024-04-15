clc; clear;
close all;

%% parameters setting %%
SorNum = 1;
c = 343;
Fs = 48000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 16000;           % 欲 resample 成的取樣頻率
MicNum = 6;           % 實驗麥克風數量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ULA %
MicStart = [0, 0, 0];
spacing = 0.07;
MicPos = zeros(MicNum, 3);
for i = 1:MicNum
    MicPos(i, :) = [MicStart(1, 1) + (i-1)*spacing, MicStart(1, 2), MicStart(1, 3)];
end

%% stft parameter %%
NFFT = 1024;
hopsize = 256;
win = hamming(NFFT);
frequency = NFFT/2 + 1;
freqs_vector = linspace(0, fs/2, frequency);

%% read mic 音檔再做 stft (Y) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Second = 28;             % 使用時間長度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorLen =  Second*Fs;
sorLen =  Second*fs;

% load y %
y = zeros(MicNum, SorLen);
for i = 1:MicNum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_filename_str = ['wav_exp\', string(i), '.wav'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_filename = join(y_filename_str, '');
    [y(i, :), ~] = audioread(y_filename, [1, SorLen]);
end

% resample %
y = resample(y, 1, Fs/fs, Dimension=2);

% y do stft get Y %
y_transpose = y.';
[Y, ~, ~] = stft(y_transpose, fs, Window=win, OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
NumOfFrame = size(Y, 2);

%% compute Ryy %%
Ryy = zeros(MicNum, MicNum, frequency);
for n = 1:frequency
    for FrameNo = 1:NumOfFrame
        Ryy(:, :, n) = Ryy(:, :, n) + squeeze(Y(n, FrameNo, :))*squeeze(Y(n, FrameNo, :))';
    end
    
end

Ryy = Ryy/NumOfFrame;

%% check angle function is convex or not %%
angle = -90:1:90;
output_abs = zeros(size(angle, 2), 1);
angle_cali = 1 - angle(:, 1);
dia_load_beamformer = 10^(-2);

for ang = angle
    kappa = [sind(ang), cosd(ang), 0];
    for n = 1:frequency
        omega = 2*pi*freqs_vector(n);
        steer_vec = exp(1j*omega/c*kappa*MicPos.').';
        w_MPDR = inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*steer_vec/(steer_vec'*inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*steer_vec);
        for FrameNo = 1:NumOfFrame
            output_abs(ang+angle_cali, :) = output_abs(ang+angle_cali, :) + abs(w_MPDR'*squeeze(Y(n, FrameNo, :)));
        end

    end

end

figure(1)
plot(output_abs.');
title('angle function')
xlabel('angle')
ylabel('output power')
shg

%% GSS for anglewise freefield plane wave localization %%
[~, max_integer_angle] = max(output_abs);
gss_left = (max_integer_angle - 1) - angle_cali;
gss_right = (max_integer_angle + 1) - angle_cali;

%% check distance function is convex or not %%



%% GSS for distancewise freefield spherical wave localization %%




