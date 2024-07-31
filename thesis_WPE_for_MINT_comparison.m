clc; clear;
close all;

fs = 16000;
look_mic = 7;

for t60 = 2:16
    y_wpe_filename_str = ['y\y_wpe_', string(t60/10), '.mat'];
    y_wpe_filename = join(y_wpe_filename_str, '');
    y_wpe = cell2mat(struct2cell(load(y_wpe_filename)));

    point_start_save = 18*fs;

    % ratio_y_wpe = 0.8 / max(abs(y_wpe(look_mic, point_start_save:end)));
    y_filemane_str = ['wav\y_wpe_partial_', string(t60/10), '.wav'];
    y_filemane = join(y_filemane_str, '');
    audiowrite(y_filemane, y_wpe(look_mic, point_start_save:end), fs)
end