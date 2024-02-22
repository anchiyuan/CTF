clc; clear; 
close all;
tic

%% RIR parameter %%
SorNum = 1;                                              % source number
MicNum = 30;                                             % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)
Ts = 1/fs;                                               % Sample period (s)

% UCA %
% radius_UCA = 1;                                        % aperture (m) 
% angle_step = 360/MicNum;
% MicPos = zeros(MicNum, 3);
% for i = 1:MicNum
%     MicPos(i, :) = [1+radius_UCA*cosd((angle_step)*(i - 1)) 3+radius_UCA*sind((angle_step)*(i - 1)) 1];
% end

% ULA %
MicStart = [1, 1.5, 1];
spacing = 0.02;
MicPos = zeros(MicNum, 3);
for i = 1 : MicNum
    MicPos(i, :) = [MicStart(1, 1)+(i-1)*spacing MicStart(1, 2) MicStart(1, 3)];
end

SorPos = [2, 2.6, 1];
room_dim = [5, 6, 2.5];                                  % Room dimensions [x y z] (m)
reverberation_time = 0.2;                                % Reverberation time (s)
points_rir = 4096;                                       % Number of rir points (需比 reverberation time 還長)
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 0;                                           % Disable high-pass filter

% 畫空間圖 %
figure(1);
plot3( [0 room_dim(1, 1) room_dim(1, 1) 0 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1)], ...
       [0 0 room_dim(1, 2) room_dim(1, 2) 0 0 0 room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) 0 0 0], ...
       [0 0 0 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0] , 'k')
hold on
plot3(MicPos(:, 1), MicPos(:, 2), MicPos(:, 3), 'r.', 'MarkerSize', 1)
hold on
plot3(SorPos(:, 1), SorPos(:, 2), SorPos(:, 3), '*', 'MarkerSize', 20)
hold off
xlabel('x\_axis')
ylabel('y\_axis')
zlabel('z\_axis')
title('空間圖')
shg

%% generate ground-truth RIR (h) %%
% 產生 RIR 和存.mat 檔 %
% h = rir_generator(c, fs, MicPos, SorPos, room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
% rir_filename_str = ['h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
% rir_filemane = join(rir_filename_str, '');
% save(rir_filemane, 'h')

% load .mat 檔 %
rir_filename_str = ['h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

% 畫 ground-truth RIR time plot %
look_mic = 10;
figure(2)
plot(h(look_mic, :), 'r');
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h')
xlabel('points')
ylabel('magnitude')
shg

%% compute ground-truth CTF (H) %%
NWIN_H = 2048;
NFFT_H = 2048;
hopsize_H = floor(NWIN_H/2);    % OverlapLength = NWIN - hopsize

% RIR 轉頻域 %
h_transpose = h(:, 1:NWIN_H).';
[H, freqs_vector, H_t_vector] = stft(h_transpose, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, FrequencyRange='onesided');
frequency = size(freqs_vector, 1);
H_dB = mag2db(abs(H));

% 畫 ground-truth RIR frequency plot %
figure(3);
semilogx(freqs_vector, H_dB(:, :, look_mic), 'r')
xlim([0 fs/2])
H_yaxis_underlimit = min(H_dB(2:end, :, look_mic), [], 'all') - 5;
H_yaxis_upperlimit = max(H_dB(2:end, :, look_mic), [], 'all') + 5;
ylim([H_yaxis_underlimit H_yaxis_upperlimit])
title('H')
xlabel('frequency(Hz)')
ylabel('magnitude(dB)')
shg

%% 讀音檔 or 產生 white noise source (source) %%
Second = 20;
SorLen =  Second*fs;

% load source %
% [source, pathname] = read_audio(SorNum, Second, fs);     % speech source
[source_transpose, fs] = audioread('245.wav', [1, SorLen]);
source = source_transpose.';
% source = wgn(1, SorLen, 0);                              % white noise source

NWIN = 2048;
NFFT = 2048;
hopsize = 64;
[S, freqs_vector, S_t_vector] = stft(source, fs, Window=hamming(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
mean_S = mean(abs(S), 'all');

%% RIR mix source 先在時域上 convolution 再做 stft (y and Y) %%
% convolution source and RIR %
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(h(i, :), source);
end

as = as(:, 1:SorLen);

% 加上 white noise 當作 interferer %
% SNR = 40;
% for i = 1:MicNum
%     y(i, :) = awgn(as(i, :), SNR, 'measured');
% end

% 不加 white noise %
y = as;

% y 轉頻域 %
y_transpose = y.';
[Y, freqs_vector, Y_t_vector] = stft(y_transpose, fs, Window=hamming(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
Y_dB = mag2db(abs(Y));
NumOfFrame = size(Y_t_vector, 1);
NumOfFrame_vector = 1:1:NumOfFrame;

%% initial a %%
% use freefield point source model initialize a %
a = zeros(MicNum, SorNum, frequency);
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((SorPos - MicPos(i, :)).^2));
end

for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, :, n) = exp(-1j*omega/c*distance)./distance;
end

% adjust a_ini_tdomain 之最高點與 h 之最高點等高 讓頻域圖高度相同 %
a_ini_sym = [squeeze(a) conj(flip(squeeze(a(:, :, 2:end-1)), 2))]; 
a_ini_tdomain = ifft(a_ini_sym, NFFT_H, 2, 'symmetric');
h_max = max(h(look_mic, :));
a_ini_tdomain_max = max(a_ini_tdomain(look_mic, :));
ratio_a_ini_tdomain = h_max/a_ini_tdomain_max;


for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, :, n) = ratio_a_ini_tdomain*exp(-1j*omega/c*distance)./distance;
end

a_ini_sym = [squeeze(a) conj(flip(squeeze(a(:, :, 2:end-1)), 2))]; 
a_ini_tdomain = ifft(a_ini_sym, NFFT_H, 2, 'symmetric');

% 畫 adjusted a_ini time plot %
figure(4);
plot(h(look_mic, :), 'r')
hold on
plot(a_ini_tdomain(look_mic, :), 'b')
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('a\_ini\_tdomain')
xlabel('points')
ylabel('magnitude')
legend('h', 'a\_ini\_tdomai')
shg

%% initial Ryy %%
ini_frame = floor(NumOfFrame/10);
gamma = 0.99;
Ryy = zeros(MicNum, MicNum, frequency);
for FrameNo =1:ini_frame
    Y_choose = squeeze(Y(:, FrameNo, :)).';
    for n = 1:frequency
         Ryy(:, :, n) = gamma*Ryy(:, :, n) + (1 - gamma)*Y_choose(:, n) * Y_choose(:, n)';  
    end

end

% Ryy = Ryy / ini_frame;

%% initial Gsy Gss %%
Gsy = zeros(MicNum, 1, frequency);
Gss = zeros(1, frequency);

%% recursive process %%
alpha = 0.99;
beta = 0.99;
beamformer_mode = 'DAS';    % 'DAS' or 'MPDR'
a_plot_mode = 'time';    % 'frequency' or 'time'
w = zeros(MicNum, 1, frequency);
s_hat = zeros(1, frequency);
shat_thres_speech = mean_S + 0.1;
shat_thres_whitenoise = 39;
save_a = zeros(NumOfFrame, MicNum, 1, frequency);
save_shat = zeros(frequency, NumOfFrame);

for FrameNo = ini_frame+1:NumOfFrame
    for n = 1:frequency
        Y_choose = squeeze(Y(n, FrameNo, :));

        %  設置 diagonal loading 來判斷使用哪種 beamformer %
        if strcmp(beamformer_mode, 'DAS')
            % update with DAS %
            dia_load_beamformer = 10000;
        elseif strcmp(beamformer_mode, 'MPDR')
            % update with MPDR %
            dia_load_beamformer = 10^(-1);
        end

        w(:, 1, n) = inv(Ryy(:, :, n) + dia_load_beamformer*eye(MicNum))*a(:, :, n) / real(a(:, :, n)'*inv(Ryy(:, :, n) + dia_load_beamformer*eye(MicNum))*a(:, :, n));
        s_hat(:, n) = w(:, 1, n)'*Y_choose;

        % 儲存 beam 到的 S_hat %
        save_shat(n, FrameNo) = s_hat(:, n);
        
        % 制止 S_hat 不夠大時的錯誤更新 %
        if abs(s_hat(:, n)) > shat_thres_speech
            Gss(1, n) = alpha*Gss(1, n) + (1 - alpha)*conj(s_hat(:, n))*s_hat(:, n);
            Gsy(:, 1, n) = alpha*Gsy(:, 1, n) + (1 - alpha)*conj(s_hat(:, n))*Y_choose;
            a(:, :, n) = Gsy(:, 1, n)./Gss(1, n);
        end

        Ryy(:, :, n) = beta*Ryy(:, :, n) + (1 - beta)*Y_choose*Y_choose';
    end

    % 儲存每次 recursive 的 a %
    save_a(FrameNo, :, :, :) = a;
    
    % 一直畫 a 的變化 %
    figure(5);
    if strcmp(a_plot_mode, 'frequency')
        a_dB = mag2db(abs(a));
    
        semilogx(freqs_vector, squeeze(a_dB(look_mic, :, :)), 'b')
        hold on
        semilogx(freqs_vector, H_dB(:, :, look_mic), 'r')
        hold off
        xlim([0 fs/2])
        ylim([H_yaxis_underlimit H_yaxis_upperlimit])
        title('a', num2str(FrameNo))
        xlabel('frequency(Hz)')
        ylabel('magnitude(dB)')
        legend('Predicted RIR', 'RIR (ground truth)')
        shg
    elseif strcmp(a_plot_mode, 'time')
        a_sym = [squeeze(a) conj(flip(squeeze(a(:, :, 2:end-1)), 2))]; 
        a_tdomain = ifft(a_sym, NFFT_H, 2, 'symmetric');

        plot(h(look_mic, :), 'r')
        hold on
        plot(a_tdomain(look_mic, :), 'b')
        hold off
        ylim([h_yaxis_underlimit h_yaxis_upperlimit])
        title('a\_tdomain', num2str(FrameNo))
        xlabel('points')
        ylabel('magnitude')
        legend('h', 'a\_tdomain')
        shg

    end

    % print 出先在處理哪一個 frame %
    fprintf('right now processing frame = %d\n', FrameNo)
end

%% 畫 a_fin frequency and time plot %%
plot_mode = 'semilogx';    % 'semilogx' or 'linear'
a_dB = mag2db(abs(a));

figure(6);
if strcmp(plot_mode, 'semilogx')
    semilogx(freqs_vector, squeeze(a_dB(look_mic, :, :)), 'b')
    hold on
    semilogx(freqs_vector, H_dB(:, :, look_mic), 'r')
    hold off
elseif strcmp(plot_mode, 'linear')
    plot(freqs_vector, squeeze(a_dB(look_mic, :, :)), 'b')
    hold on
    plot(freqs_vector, H_dB(:, :, look_mic), 'r')
    hold off
end

xlim([0 fs/2])
ylim([H_yaxis_underlimit H_yaxis_upperlimit])
title('a\_fin')
xlabel('frequency(Hz)')
ylabel('magnitude(dB)')
legend('Predicted RIR', 'RIR (ground truth)')
shg

% 轉回時域畫圖比較 %
a_fin_sym = [squeeze(a) conj(flip(squeeze(a(:, :, 2:end-1)), 2))]; 
a_fin_tdomain = ifft(a_fin_sym, NFFT_H, 2, 'symmetric');

h_max = max(h(look_mic, :));
a_fin_tdomain_max = max(a_fin_tdomain(look_mic, :));
ratio_a_fin_tdomain = h_max/a_fin_tdomain_max;


figure(7);
plot(h(look_mic, :), 'r')
hold on
plot(a_fin_tdomain(look_mic, :)*ratio_a_fin_tdomain, 'b')
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('a\_fin\_tdomain')
xlabel('points')
ylabel('magnitude')
legend('h', 'a\_fin\_tdomain')
shg

%% plot spectrogram (S save_shat) %%
S_dB = mag2db(abs(S));

figure(8);
mesh(NumOfFrame_vector, freqs_vector, S_dB)
S_caxis_upperlimit = max(S_dB, [], 'all') + 5;
S_caxis_underlimit = -90;
colorbar
clim([S_caxis_underlimit, S_caxis_upperlimit])
view(2)
xlim([1 NumOfFrame])
title('S')
xlabel('frame')
ylabel('frequency(Hz)')
shg

save_shat_dB = mag2db(abs(save_shat));

figure(9);
mesh(NumOfFrame_vector, freqs_vector, save_shat_dB)
colorbar
clim([S_caxis_underlimit, S_caxis_upperlimit])
view(2)
xlim([1 NumOfFrame])
title('save\_shat')
xlabel('frame')
ylabel('frequency(Hz)')

%% save_shat 還原回時域 (shat_tdomain) %%
% save_shat 轉回時域 %
[shat_tdomain_transpose] = istft(save_shat, fs, Window=hamming(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, ConjugateSymmetric=true, FrequencyRange='onesided');
shat_tdomain = shat_tdomain_transpose.';

% adjust shat_tdomain 的最高點使之與 source 的一樣 %
source_max  = max(abs(source(1, (ini_frame*hopsize+NWIN):end)));
shat_tdomain_max  = max(abs(shat_tdomain));
ratio_shat_tdomain = source_max/shat_tdomain_max;
shat_tdomain = shat_tdomain.*ratio_shat_tdomain;

% 畫 shat_tdomain time plot %
figure(10)
plot(source(1, :), 'r');
hold on
plot(shat_tdomain(1, :), 'b');
hold off
title('shat\_tdomain')
xlabel('points')
ylabel('magnitude')
legend('source', 'shat\_tdomain')
shg

% save .wav 檔 %
audiowrite('source.wav', source(1, :), fs)
ratio_y = 1 / max(abs(y(look_mic, :))) ;
audiowrite('y.wav', y(look_mic, :)*ratio_y, fs)
shat_tdomain_filename_str = ['shat_tdomain-beamformer_mode=', string(beamformer_mode), '-alpha=', string(alpha), '.wav'];
shat_tdomain_filemane = join(shat_tdomain_filename_str, '');
audiowrite(shat_tdomain_filemane, shat_tdomain, fs)

%% save .mat 檔 %%
% a_filename_str = ['a_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
% a_filemane = join(a_filename_str, '');
% save(a_filemane,'a')
