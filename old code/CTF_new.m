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
% radius_UCA = 0.5;                                        % aperture (m) 
% angle_step = 360/MicNum;
% MicPos = zeros(MicNum, 3);
% for mic_num = 1 : MicNum
%     MicPos(mic_num, :) = [1 + radius_UCA*cosd((angle_step)*(mic_num - 1)) 3 + radius_UCA*sind((angle_step)*(mic_num - 1)) 1];
% end

% ULA %
MicStart = [1, 1.5, 1];
spacing = 0.02;
MicPos = zeros(MicNum, 3);
for i = 1 : MicNum
    MicPos(i, :) = [MicStart(1,1)+(i-1)*spacing MicStart(1,2) MicStart(1,3)];
end

SorPos = [2, 2.6, 1];
room_dim = [5, 6, 2.5];                                  % Room dimensions [x y z] (m)
reverberation_time = 0.4;                                % Reverberation time (s)
points_rir = 8192;                                       % Number of rir points (需比 reverberation time 還長)
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 0;                                           % Disable high-pass filter

% 畫空間圖 %
% figure(1);
% plot3( [0 room_dim(1,1) room_dim(1,1) 0 0 0 room_dim(1,1) room_dim(1,1) 0 0 room_dim(1,1) room_dim(1,1) 0 0 room_dim(1,1) room_dim(1,1)], ...
%        [0 0 room_dim(1,2) room_dim(1,2) 0 0 0 room_dim(1,2) room_dim(1,2) room_dim(1,2) room_dim(1,2) room_dim(1,2) room_dim(1,2) 0 0 0], ...
%        [0 0 0 0 0 room_dim(1,3) room_dim(1,3) room_dim(1,3) room_dim(1,3) 0 0 room_dim(1,3) room_dim(1,3) room_dim(1,3) room_dim(1,3) 0] , 'k')
% hold on
% plot3(MicPos(:,1), MicPos(:,2), MicPos(:,3), 'r.', 'MarkerSize', 1)
% hold on
% plot3(SorPos(:,1), SorPos(:,2), SorPos(:,3), '*', 'MarkerSize', 20)
% hold off
% xlabel('x\_axis')
% ylabel('y\_axis')
% zlabel('z\_axis')
% title('空間圖')
% shg

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
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h')
xlabel('points')
ylabel('magnitude')
shg

%% compute ground-truth CTF (H) %%
NWIN_H = 512;
NFFT_H = 1024;
hopsize_H = 128;    % OverlapLength = NWIN - hopsize

% RIR 轉頻域 %
h_transpose = h.';
[H, freqs_vector, H_t_vector] = stft(h_transpose, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, FrequencyRange='onesided');
frequency = size(freqs_vector, 1);
L = size(H_t_vector, 1);
L_vector = 1:1:L;
H_dB = mag2db(abs(H));

% 畫 ground-truth RIR frequency plot %
figure(3);
mesh(L_vector, freqs_vector, H_dB(:, :, look_mic))
H_caxis_upperlimit = max(H_dB(:, :, look_mic), [], 'all') + 5;
H_caxis_underlimit = min(H_dB(:, :, look_mic), [], 'all') - 5;
colorbar
clim([H_caxis_underlimit, H_caxis_upperlimit])
view(2)
title('H')
xlabel('frame')
ylabel('frequency(Hz)')
shg

% figure(23);
% plot(freqs_vector, squeeze(H_dB(:, :, look_mic)).')
% ylim([H_caxis_underlimit, H_caxis_upperlimit])
% title('H')
% xlabel('frequency(Hz)')
% ylabel('dB')
% shg

%% reconstruct RIR from H (h_recon) %%
% H 轉回時域 %
% [h_recon_transpose] = istft(H, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, FrequencyRange='onesided');
% h_recon = h_recon_transpose.';

% 畫 reconstruct RIR time plot %
% figure(4)
% plot(h(look_mic, :), 'r');
% hold on
% plot(h_recon(look_mic, :), 'b');
% hold off
% ylim([h_yaxis_underlimit h_yaxis_upperlimit])
% title('h\_recon')
% xlabel('points')
% ylabel('magnitude')
% legend('h', 'h\_recon')
% shg

%% 讀音檔 or 產生 white noise source (source) %%
Second = 20;
SorLen =  Second*fs;

% load source %
% [source, pathname] = read_audio(SorNum, Second, fs);     % speech source
[source_transpose, fs] = audioread('245.wav', [1, SorLen]);
source = source_transpose.';
% source = wgn(1, SorLen, 0);                              % white noise source

%% compute source signal for frequency (S) %%
NWIN = 1024;
NFFT = NFFT_H;
hopsize = hopsize_H;    % OverlapLength = NWIN - hopsize

% source 轉頻域 %
source_transpose = source.';
[S, freqs_vector, S_t_vector] = stft(source_transpose, fs, Window=hamming(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
NumOfFrame = size(S_t_vector, 1);
NumOfFrame_vector = 1:1:NumOfFrame;
S_dB = mag2db(abs(S));

% 畫 source frequency plot %
figure(5);
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

% 利用 CTF model 還原 mic frequency stft (Y_CTF) %
Y_CTF = zeros(frequency, NumOfFrame, MicNum);
for FrameNo= L:NumOfFrame
    for i = 1 : MicNum
        Y_CTF(:, FrameNo, i) = sum(H(:, :, i).*flip(S(:, FrameNo-L+1:FrameNo), 2), 2);
    end
end

% Y_CTF 轉時域 (y_ctf) %
[y_CTF_transpose] = istft(Y_CTF, fs, Window=rectwin(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, ConjugateSymmetric=true, FrequencyRange='onesided');
y_CTF = y_CTF_transpose.';


% 畫 mic signal 比較 time plot %
ratio_y = 1 / max(abs(y(look_mic, :))) ;
ratio_y_CTF = 1 / max(abs(y_CTF(look_mic, :))) ;

figure(7);
plot(y(look_mic, :)*ratio_y, 'r');
hold on
plot(y_CTF(look_mic, :)*ratio_y_CTF, 'b');
hold off
ylim([-1.1 1.1])
title('y and y\_ctf')
xlabel('points')
ylabel('magnitude')
legend('y', 'y\_ctf')
shg

% 存兩種 .wav 檔 %
audiowrite('wav\y.wav', y(look_mic, :)*ratio_y, fs)
audiowrite('wav\y_CTF.wav', y_CTF(look_mic, :)*ratio_y_CTF, fs)

%% initial Rss Rsy %%
ini_frame = floor(NumOfFrame/3);

% 初始化 Rss %
forfac_ini_Rss = 0.999;
Rss = zeros(L, L, frequency);
for FrameNo= L:ini_frame
    S_choose = S(:, FrameNo+1-L:FrameNo).';
    for n = 1: frequency
         Rss(:,:,n) = forfac_ini_Rss*Rss(:,:,n) + (1-forfac_ini_Rss)*flip(S_choose(:,n)) * flip(S_choose(:,n))';
    end  

end

% 初始化 Rsy %
forfac_ini_Rsy = 0.999;
Rsy = zeros(L, MicNum, frequency);
for FrameNo= L:ini_frame
    S_choose = S(:, FrameNo+1-L:FrameNo).';
    Y_choose = squeeze(Y(:, FrameNo, :)).';
    for n = 1: frequency
         Rsy(:,:,n) = forfac_ini_Rsy*Rsy(:,:,n) + (1-forfac_ini_Rsy)*flip(S_choose(:,n)) * Y_choose(:,n)';
    end 

end

%% 畫圖看初始 A 在頻域的樣子 (A_ini) %%
A_ini = zeros(MicNum, L, frequency);
dia_load_ini = 10^(-7);
for n = 1:frequency
    A_ini(:,:,n) = Rsy(:,:,n)'*inv(Rss(:,:,n) + dia_load_ini.*eye(L));
end

A_ini_dB = mag2db(abs(A_ini));

% 畫 A_ini frequency plot %
figure(8);
mesh(L_vector, freqs_vector, squeeze(A_ini_dB(look_mic, :, :)).')
colorbar
clim([H_caxis_underlimit, H_caxis_upperlimit])
view(2)
xlim([1 L])
title('A\_ini')
xlabel('frame')
ylabel('frequency(Hz)')
shg

% figure(28);
% plot(freqs_vector, squeeze(A_ini_dB(look_mic, :, :)))
% ylim([H_caxis_underlimit, H_caxis_upperlimit])
% title('A\_ini')
% xlabel('frequency(Hz)')
% ylabel('dB')
% shg

%% 畫圖看初始 A 在時域的樣子 有adjust最高點 (A_ini_tdomain) %%
% A_ini 轉回時域 %
A_ini_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_ini_forplot(:, :, i) = squeeze(A_ini(i, :, :)).';
end

[A_ini_tdomain_transpose] = istft(A_ini_forplot, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, ConjugateSymmetric=true, FrequencyRange='onesided');
A_ini_tdomain = A_ini_tdomain_transpose.';

% adjust 圖中的 A_ini_tdomain 使最高點與 h 一樣 %
look_mic = 10;    % 在宣告一次 用來看其他mic
h_max = max(h(look_mic, :));
A_ini_tdomain_max = max(A_ini_tdomain(look_mic, :));
ratio_A_ini_tdomain = h_max/A_ini_tdomain_max;

% 畫 A_ini time plot %
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;

figure(9)
plot(h(look_mic, :), 'r');
hold on
% plot(A_ini_tdomain(look_mic, :).*ratio_A_ini_tdomain, 'b');
plot(A_ini_tdomain(look_mic, :), 'b');
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('A\_ini\_tdomain')
xlabel('points')
ylabel('magnitude')
legend('h', 'A\_ini\_tdomain')
shg

% 算 h 與 A_ini_tdomain 的 matching error %
single_proof = 10^(-6);
rms = sum((h(look_mic, :) - A_ini_tdomain(look_mic, :)).^2)/sum((h(look_mic, :) + single_proof).^2);

cosine_similarity = dot(h(look_mic, :), A_ini_tdomain(look_mic, :)) / (norm(h(look_mic, :), 2)*norm(A_ini_tdomain(look_mic, :), 2));

%% initial A RAA rAy %%
A = zeros(MicNum, L, frequency);
RAA = zeros(L, L, frequency);
rAy = zeros(L, 1, frequency);

%% recursive process %%
alpha = 0.999;
beta = 0.999;
gamma = 0.001;
delta = 0.001;
dia_load_A = 10^(-7);
dia_load_S_predict = 10^(-8);
RssRsy_mode = 'number';    % 'number' or 'framecountbased'
A_mode = 'true';    % 'predict' or 'true'
save_mode = 'front';    % 'front' or 'back'
S_predict = zeros(L, 1, frequency);
S_save = zeros(frequency, NumOfFrame);
S_save_middle = zeros(frequency, NumOfFrame);
S_save_L = zeros(frequency, NumOfFrame);
A_forplot = zeros(frequency, L, MicNum);

process_first_frame = ini_frame+1;
for FrameNo = process_first_frame : NumOfFrame
    for n = 1:frequency
        Y_before = squeeze(Y(n, FrameNo-1, :));
        Y_now = squeeze(Y(n, FrameNo, :));
        if strcmp(A_mode, 'predict')
            if FrameNo ~= process_first_frame
                if strcmp(RssRsy_mode, 'number')
                    % update Rss Rsy using alpha and beta %
                    Rss(:,:,n) = alpha.*Rss(:,:,n) + (1 - alpha).*S_predict(:,:,n)*S_predict(:,:,n)';
                    Rsy(:,:,n) = beta.*Rsy(:,:,n) + (1 - beta).*S_predict(:,:,n)*Y_before';
                elseif strcmp(RssRsy_mode, 'framecountbased')
                    % update Rss Rsy using frame count based %
                    Rss(:, :, n) = ((FrameNo-L-1).*Rss(:, :, n) + S_predict(:, :, n)*S_predict(:, :, n)')./(FrameNo-L);
                    Rsy(:, :, n) = ((FrameNo-L-1).*Rsy(:, :, n) + S_predict(:, :, n)*Y_before')./(FrameNo-L);
                end
    
            end

        end

        if strcmp(A_mode, 'predict')
            % 使用估到的 S 來更新 A %
            A(:, : ,n) = Rsy(:, :, n)'*inv(Rss(:, :, n) + dia_load_A.*eye(L));
        elseif strcmp(A_mode, 'true')
            % 使用 ground-truth H 來當作 A %
            for i = 1 : MicNum
                A(i, :, :) = squeeze(H(:, :, i)).';
            end

        end

        RAA(:, :, n) = gamma.*RAA(:, :, n) + (1 - gamma).*A(:, :, n)'*A(:, :, n);
        rAy(:, :, n) = delta.*rAy(:, :, n) + (1 - delta).*A(:, :, n)'*Y_now;
        S_predict(:, :, n) = inv(RAA(:, :, n) + dia_load_S_predict.*eye(L))*rAy(:, :, n);
    end
    
    % 儲存 S_predict %
    if strcmp(save_mode, 'front')
        % 存最前面的 frame %
        S_save(:, FrameNo) = S_predict(1, :, :);
        S_save_middle(:, FrameNo) = S_predict(floor(L/2), :, :);
        S_save_L(:, FrameNo) = S_predict(L, :, :);
    elseif strcmp(save_mode, 'back')
        % 存最後面的 frame %
        S_save(:, FrameNo-(L-1)) = S_predict(L, :, :);
        if FrameNo == NumOfFrame
            for count = 1:L
                S_save(:, NumOfFrame+1-count) = S_predict(count, :, :);
            end
    
        end

    end
    
    % 一直畫 A_tdomain 的變化 %
    if strcmp(A_mode, 'predict')
        for i = 1 : MicNum
            A_forplot(:, :, i) = squeeze(A(i, :, :)).';
        end
        
        [A_tdomain_transpose] = istft(A_forplot, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, ConjugateSymmetric=true, FrequencyRange='onesided');
        A_tdomain = A_tdomain_transpose.';
         
        figure(13)
        plot(h(look_mic, :), 'r');
        hold on
        plot(A_tdomain(look_mic, :), 'b');
        hold off
        ylim([h_yaxis_underlimit h_yaxis_upperlimit])
        title('A\_tdomain', num2str(FrameNo))
        xlabel('points')
        ylabel('magnitude')
        legend('h', 'A\_tdomain')
        shg
    end

    % print 出先在處理哪一個 frame %
    fprintf('right now processing frame = %d\n', FrameNo)
end

%% 畫圖看最後 A 在時域的樣子 有adjust最高點 (A_fin_tdomain) %%
% A_fin 轉回時域 %
A_fin_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_fin_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

[A_fin_tdomain_transpose] = istft(A_fin_forplot, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, ConjugateSymmetric=true, FrequencyRange='onesided');
A_fin_tdomain = A_fin_tdomain_transpose.';

% adjust 圖中的 A_fin_tdomain 使最高點與 h 一樣 %
look_mic = 10;    % 在宣告一次 用來看其他mic
h_max = max(h(look_mic, :));
A_fin_tdomain_max = max(A_fin_tdomain(look_mic, :));
ratio_A_fin_tdomain = h_max/A_fin_tdomain_max;

% 畫 A_fin time plot %
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;

figure(14)
plot(h(look_mic, :), 'r');
hold on
plot(A_fin_tdomain(look_mic, :), 'b');
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('A\_fin\_tdomain')
xlabel('points')
ylabel('magnitude')
legend('h', 'A\_ini\_tdomain')
shg

%% plot spectrogram (S_save S_save_L) %%
S_save_dB = mag2db(abs(S_save));

figure(15);
mesh(NumOfFrame_vector, freqs_vector, S_save_dB)
colorbar
clim([S_caxis_underlimit, S_caxis_upperlimit])
view(2)
xlim([1 NumOfFrame])
title('S\_save')
xlabel('frame')
ylabel('frequency(Hz)')

S_save_L_dB = mag2db(abs(S_save_L));

figure(16);
mesh(NumOfFrame_vector, freqs_vector, S_save_L_dB)
colorbar
clim([S_caxis_underlimit, S_caxis_upperlimit])
view(2)
xlim([1 NumOfFrame])
title('S\_save\_L')
xlabel('frame')
ylabel('frequency(Hz)')

%% predicted source 還原回時域 (source_predict) %%
% S_save 轉回時域 %
[source_predict_transpose] = istft(S_save, fs, Window=hamming(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, ConjugateSymmetric=true, FrequencyRange='onesided');
source_predict = source_predict_transpose.';

% adjust source_predict 的最高點使之與 source 的一樣 %
source_max  = max(abs(source(1, (ini_frame*hopsize+NWIN):end)));
source_predict_max  = max(abs(source_predict(1, (ini_frame*hopsize+NWIN):(Second-0.1)*fs)));
ratio_source_predict = source_max/source_predict_max;
source_predict = source_predict.*ratio_source_predict;

% 畫 source_predict time plot %
figure(17)
plot(source(1, :), 'r');
hold on
plot(source_predict(1, :), 'b');
hold off
ylim([source_yaxis_underlimit source_yaxis_upperlimit])
title('source\_predict')
xlabel('points')
ylabel('magnitude')
legend('source', 'source\_predict')
shg

%% save .wav 檔 %%
audiowrite('wav\source.wav', source(1, (ini_frame*hopsize+NWIN):end), fs)
ratio_y = 1 / max(abs(y(look_mic, (ini_frame*hopsize+NWIN):end))) ;
audiowrite('wav\y.wav', y(look_mic, (ini_frame*hopsize+NWIN):end)*ratio_y, fs)
source_predict_filemane_str = ['wav\source_predict-', string(reverberation_time), 'x', string(NWIN_H), 'x', string(NFFT_H), 'x', string(hopsize_H), 'x',...
    string(NWIN), 'x', string(NFFT), 'x',string(hopsize), '_A_mode=', A_mode, '.wav'];
source_predict_filemane = join(source_predict_filemane_str, '');
audiowrite(source_predict_filemane, source_predict(1, (ini_frame*hopsize+NWIN):end), fs)
fprintf(source_predict_filemane)
fprintf('\n')

toc
