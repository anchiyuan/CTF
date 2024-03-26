clc; clear;
close all;

fs = 48000;
second = 30;
source = wgn(1, fs*second, 1, 'linear');    % wgn(row, column, variance, 'linear');
max_mag = max(abs(source));
source = source/max_mag;
filename_str = ['white_noise_', string(second), 's.wav'];
filename = join(filename_str, '');
audiowrite(filename, source, fs)

NFFT = 2048;
hopsize = NFFT/4;
[S, ~, ~] = stft(source.', fs, Window=hann(NFFT), OverlapLength=NFFT-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
S_dB = mag2db(abs(S));
figure(3)
mesh(S_dB)
colorbar
view(2)