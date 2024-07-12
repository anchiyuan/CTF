clc; clear;
close all;

fs = 16000;

algorithm = ["Wiener", "RLS", "Kalman"];
reverberation_time = 1.6;
MicNum = 38;
points_rir = 32768;
look_mic = 38;

rir_filename_str = ['h\h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

f1 = figure;
f2 = figure;

for i = 1:size(algorithm, 2)
    al = algorithm(:, i);
    A_tdomain_filename_str = ['A_tdomain\CTF_combined_', string(reverberation_time), 'x', al, '_A_tdomain.mat'];
    A_tdomain_filename = join(A_tdomain_filename_str, '');
    A_tdomain = cell2mat(struct2cell(load(A_tdomain_filename)));
    
    ATF = fft(h, points_rir, 2);
    ATF_estimated = fft(A_tdomain, points_rir, 2);
    
    % 畫 A_tdomain time plot %
    figure(f1)
    subplot(size(algorithm, 2), 1, i);
    plot(h(look_mic, :), 'r');
    hold on
    plot(A_tdomain(look_mic, :), 'b-.');
    hold off
    h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
    h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
    ylim([h_yaxis_underlimit h_yaxis_upperlimit])
    title(al)
    legend('ground-truth RIR', 'estimated RIR')
    xlabel('points')
    ylabel('amplitude')
    
    figure(f2)
    subplot(size(algorithm, 2), 2, i*2-1);
    semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF(look_mic, 1:points_rir/2+1))), 'r');
    hold on
    semilogx(linspace(0, fs/2, points_rir/2+1), 20*log10(abs(ATF_estimated(look_mic, 1:points_rir/2+1))), 'b-.');
    hold off
    xlim([200 8000])
    title(al)
    legend('ground-truth ATF', 'estimated ATF')
    xlabel('frequency (Hz)')
    ylabel('dB')
    
    subplot(size(algorithm, 2), 2, i*2);
    semilogx(linspace(0, fs/2, points_rir/2+1), unwrap(angle(ATF(look_mic, 1:points_rir/2+1))), 'r');
    hold on
    semilogx(linspace(0, fs/2, points_rir/2+1), unwrap(angle(ATF_estimated(look_mic, 1:points_rir/2+1))), 'b-.');
    hold off
    xlim([200 8000])
    title(al)
    legend('ground-truth ATF', 'estimated ATF')
    xlabel('frequency (Hz)')
    ylabel('phase (radius)')

end

