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
reverberation_time = 0.6;                                % Reverberation time (s)
points_rir = 14336;                                       % Number of rir points (需比 reverberation time 還長)
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
NWIN_H = 64;
NFFT_H = 2048;
hopsize_H = 32;    % OverlapLength = NWIN - hopsize

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

% figure(4);
% plot(freqs_vector, squeeze(H_dB(:, :, look_mic)).')
% ylim([H_caxis_underlimit, H_caxis_upperlimit])
% title('H')
% xlabel('frequency(Hz)')
% ylabel('dB')
% shg

%% reconstruct RIR from H (h_recon) %%
% H 轉回時域 %
% [h_recon_transpose] = istft(H, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, ConjugateSymmetric=true, FrequencyRange='onesided');
% h_recon = h_recon_transpose.';

% 畫 reconstruct RIR time plot %
% figure(5)
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
% [source, pathname] = read_audio(SorNum, Second, fs);     % speech source load with .m檔

[source_transpose, fs] = audioread('245.wav', [1, SorLen]);    % speech source load with audioread
source = source_transpose.';

% source = wgn(1, SorLen, 0);                              % white noise source

%% compute source signal for frequency (S) %%
NWIN = 64;
NFFT = NFFT_H;
hopsize = hopsize_H;    % OverlapLength = NWIN - hopsize

% source 轉頻域 %
source_transpose = source.';
[S, freqs_vector, S_t_vector] = stft(source_transpose, fs, Window=hamming(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
NumOfFrame = size(S_t_vector, 1);
NumOfFrame_vector = 1:1:NumOfFrame;
S_dB = mag2db(abs(S));

% 畫 source frequency plot %
figure(6);
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

% 畫 source and mic time plot %
% figure(7)
% plot(source(1, :), 'r');
% hold on
% plot(y(look_mic, :), 'b');
% hold off
source_yaxis_upperlimit = max(source) + 0.1;
source_yaxis_underlimit = min(source) - 0.1;
% ylim([source_yaxis_underlimit source_yaxis_upperlimit])
% title('source and y')
% xlabel('points')
% ylabel('magnitude')
% legend('source', 'y')
% shg

% y 轉頻域 %
y_transpose = y.';
[Y, freqs_vector, Y_t_vector] = stft(y_transpose, fs, Window=hamming(NWIN), OverlapLength=NWIN-hopsize, FFTLength=NFFT, FrequencyRange='onesided');
Y_dB = mag2db(abs(Y));

% 畫 mic frequency plot %
% figure(8);
% mesh(NumOfFrame_vector, freqs_vector, Y_dB(:, :, look_mic))
% colorbar
% clim([S_caxis_underlimit, S_caxis_upperlimit])
% view(2)
% xlim([1 NumOfFrame])
% title('Y')
% xlabel('frame')
% ylabel('frequency(Hz)')
% shg

%% use recursive_ATF to generate initial S for Rss Rsy (S_hat_save) %%
RssRsy_ini_mode = 'S';    % 'S_hat_save' or 'S'

if strcmp(RssRsy_ini_mode, 'S')
    ini_frames_Ryy = floor(NumOfFrame/10);
    start_frame_Ryy = 1;
    end_frame_Ryy = start_frame_Ryy + (ini_frames_Ryy - 1);
    stable_frame_a = floor(NumOfFrame/4);
    start_frame_a = end_frame_Ryy + 1;
    ini_frames_RssRsy = floor(NumOfFrame/3);
    end_frame_a = (stable_frame_a + 1) + (ini_frames_RssRsy - 1);

end

if strcmp(RssRsy_ini_mode, 'S_hat_save')
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
    % figure(9);
    % plot(h(look_mic, :), 'r')
    % hold on
    % plot(a_ini_tdomain(look_mic, :), 'b')
    % hold off
    % ylim([h_yaxis_underlimit h_yaxis_upperlimit])
    % title('a\_ini\_tdomain')
    % xlabel('points')
    % ylabel('magnitude')
    % legend('h', 'a\_ini\_tdomai')
    % shg
    
    % initial Ryy %
    ini_frames_Ryy = floor(NumOfFrame/10);
    start_frame_Ryy = 1;
    end_frame_Ryy = start_frame_Ryy + (ini_frames_Ryy - 1);
    forfac_ini_Ryy = 0.99;
    Ryy = zeros(MicNum, MicNum, frequency);
    for FrameNo = start_frame_Ryy : end_frame_Ryy
        Y_choose = squeeze(Y(:, FrameNo, :)).';
        for n = 1:frequency
             Ryy(:, :, n) = forfac_ini_Ryy * Ryy(:, :, n) + (1 - forfac_ini_Ryy)*Y_choose(:, n)*Y_choose(:, n)';  
        end
    
    end
    
    % initial Gsy Gss %
    Gsy = zeros(MicNum, 1, frequency);
    Gss = zeros(1, frequency);
    
    % recursive process %
    stable_frame_a = floor(NumOfFrame/4);
    start_frame_a = end_frame_Ryy + 1;
    ini_frames_RssRsy = floor(NumOfFrame/3);
    end_frame_a = (stable_frame_a + 1) + (ini_frames_RssRsy - 1);
    beamformer_mode = 'DAS';    % 'DAS' or 'MPDR'
    forfac_GssGsy = 0.99;
    forfac_Ryy = 0.99;
    w = zeros(MicNum, 1, frequency);
    S_hat = zeros(frequency, 1);
    mean_S = mean(abs(S), 'all');
    S_hat_thres_speech = mean_S + 0.1;
    S_hat_save= zeros(frequency, NumOfFrame);
    a_plot_mode = 'time';    % 'frequency' or 'time'
    
    for FrameNo = start_frame_a : end_frame_a
        for n = 1:frequency
            Y_choose = squeeze(Y(n, FrameNo, :));
    
            %  設置 diagonal loading 來判斷使用哪種 beamformer %
            if strcmp(beamformer_mode, 'DAS')
                % update with DAS %
                dia_load_beamformer = 1000;
            elseif strcmp(beamformer_mode, 'MPDR')
                % update with MPDR %
                dia_load_beamformer = 10^(-1);
            end
    
            w(:, 1, n) = inv(Ryy(:, :, n) + dia_load_beamformer*eye(MicNum))*a(:, :, n) / real(a(:, :, n)'*inv(Ryy(:, :, n) + dia_load_beamformer*eye(MicNum))*a(:, :, n));
            S_hat(n, :) = w(:, 1, n)'*Y_choose;
            
            % 儲存 beam 到的 S_hat %
            S_hat_save(n, FrameNo) = S_hat(n, :);
    
            % 制止 S_hat 不夠大時的錯誤更新 %
            if abs(S_hat(n, :)) > S_hat_thres_speech
                Gss(1, n) = forfac_GssGsy * Gss(1, n) + (1 - forfac_GssGsy)*conj(S_hat(n, :))*S_hat(n, :);
                Gsy(:, 1, n) = forfac_GssGsy * Gsy(:, 1, n) + (1 - forfac_GssGsy)*conj(S_hat(n, :))*Y_choose;
                a(:, :, n) = Gsy(:, 1, n)./Gss(1, n);
            end
    
            Ryy(:, :, n) = forfac_Ryy * Ryy(:, :, n) + (1 - forfac_Ryy)*Y_choose*Y_choose';
        end
    
        % 一直畫 a 的變化 %
        figure(10);
        if strcmp(a_plot_mode, 'frequency')
            a_dB = mag2db(abs(a));
        
            semilogx(freqs_vector, squeeze(a_dB(look_mic, :, :)), 'b')
            hold on
            semilogx(freqs_vector, squeeze(H_dB(:, 1, look_mic)), 'r')
            hold off
            xlim([0 fs/2])
            ylim([H_caxis_underlimit H_caxis_upperlimit])
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
    
    % 畫 S_hat_save frequency plot %
    S_hat_save_dB = mag2db(abs(S_hat_save));
    
    figure(11);
    mesh(NumOfFrame_vector, freqs_vector, S_hat_save_dB)
    colorbar
    clim([S_caxis_underlimit, S_caxis_upperlimit])
    view(2)
    xlim([1 NumOfFrame])
    title('save\_shat')
    xlabel('frame')
    ylabel('frequency(Hz)')
    shg
end

%% initial Rss Rsy %%
start_frame_RssRsy = (stable_frame_a + 1) + (L - 1);
end_frame_RssRsy = end_frame_a;
% 初始化 Rss %
forfac_ini_Rss = 0.999;
Rss = zeros(L, L, frequency);
for FrameNo= start_frame_RssRsy : end_frame_RssRsy
    if strcmp(RssRsy_ini_mode, 'S_hat_save')
        S_choose = S_hat_save(:, FrameNo+1-L:FrameNo).';
    elseif strcmp(RssRsy_ini_mode, 'S')
        S_choose = S(:, FrameNo+1-L:FrameNo).';
    end

    for n = 1: frequency
         Rss(:,:,n) = forfac_ini_Rss * Rss(:,:,n) + (1 - forfac_ini_Rss)*flip(S_choose(:,n)) * flip(S_choose(:,n))';
    end

end

% 初始化 Rsy %
forfac_ini_Rsy = 0.999;
Rsy = zeros(L, MicNum, frequency);
for FrameNo= start_frame_RssRsy : end_frame_a
    if strcmp(RssRsy_ini_mode, 'S_hat_save')
        S_choose = S_hat_save(:, FrameNo+1-L:FrameNo).';
    elseif strcmp(RssRsy_ini_mode, 'S')
        S_choose = S(:, FrameNo+1-L:FrameNo).';
    end

    Y_choose = squeeze(Y(:, FrameNo, :)).';
    for n = 1: frequency
         Rsy(:,:,n) = forfac_ini_Rsy * Rsy(:,:,n) + (1 - forfac_ini_Rsy)*flip(S_choose(:,n)) * Y_choose(:,n)';  
    end

end

%% 畫圖看初始 A 在頻域的樣子 (A_ini) %%
A_ini = zeros(MicNum, L, frequency);
dia_load_A_ini = 10^(-7);
for n = 1:frequency
    A_ini(:,:,n) = Rsy(:,:,n)'*inv(Rss(:,:,n) + dia_load_A_ini.*eye(L));
end

A_ini_dB = mag2db(abs(A_ini));

% 畫 A_ini frequency plot %
figure(12);
mesh(L_vector, freqs_vector, squeeze(A_ini_dB(look_mic, :, :)).')
colorbar
clim([H_caxis_underlimit, H_caxis_upperlimit])
view(2)
xlim([1 L])
title('A\_ini')
xlabel('frame')
ylabel('frequency(Hz)')
shg

% figure(13);
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

figure(14)
plot(h(look_mic, :), 'r');
hold on
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
forfac_Rss = 0.999;
forfac_Rsy = 0.999;
forfac_RAA = 0.001;
forfac_rAy = 0.001;
dia_load_A = 10^(-7);
dia_load_S_predict = 10^(-8);
RssRsy_mode = 'forfac';    % 'forfac' or 'framecountbased'
A_mode = 'RssRsy';    % 'RssRsy' or 'H'
save_mode = 'front';    % 'front' or 'back'
S_predict = zeros(L, 1, frequency);
S_save = zeros(frequency, NumOfFrame);
S_save_L = zeros(frequency, NumOfFrame);
A_forplot = zeros(frequency, L, MicNum);

start_frame_A = end_frame_RssRsy + 1;
end_frame_A = NumOfFrame;
for FrameNo = start_frame_A : end_frame_A
    for n = 1:frequency
        Y_before = squeeze(Y(n, FrameNo-1, :));
        Y_now = squeeze(Y(n, FrameNo, :));
        if FrameNo ~= start_frame_A
            if strcmp(RssRsy_mode, 'forfac')
                % update Rss Rsy using alpha and beta %
                Rss(:,:,n) = forfac_Rss.*Rss(:,:,n) + (1 - forfac_Rss).*S_predict(:,:,n)*S_predict(:,:,n)';
                Rsy(:,:,n) = forfac_Rsy.*Rsy(:,:,n) + (1 - forfac_Rsy).*S_predict(:,:,n)*Y_before';
            elseif strcmp(RssRsy_mode, 'framecountbased')
                % update Rss Rsy using frame count based %
                Rss(:, :, n) = ((FrameNo-L-1).*Rss(:, :, n) + S_predict(:, :, n)*S_predict(:, :, n)')./(FrameNo-L);
                Rsy(:, :, n) = ((FrameNo-L-1).*Rsy(:, :, n) + S_predict(:, :, n)*Y_before')./(FrameNo-L);
            end

        end

        
        if strcmp(A_mode, 'RssRsy')
            % 使用估到的 S 來更新 A %
            A(:, : ,n) = Rsy(:, :, n)'*inv(Rss(:, :, n) + dia_load_A.*eye(L));
        elseif strcmp(A_mode, 'H')
            % 使用 ground-truth H 來當作 A %
            for i = 1 : MicNum
                A(i, :, :) = squeeze(H(:, :, i)).';
            end

        end

        RAA(:, :, n) = forfac_RAA.*RAA(:, :, n) + (1 - forfac_RAA).*A(:, :, n)'*A(:, :, n);
        rAy(:, :, n) = forfac_rAy.*rAy(:, :, n) + (1 - forfac_rAy).*A(:, :, n)'*Y_now;
        S_predict(:, :, n) = inv(RAA(:, :, n) + dia_load_S_predict.*eye(L))*rAy(:, :, n);
    end
    
    % 儲存 S_predict %
    if strcmp(save_mode, 'front')
        % 存最前面的 frame %
        S_save(:, FrameNo) = S_predict(1, :, :);
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
    if strcmp(A_mode, 'RssRsy')
        for i = 1 : MicNum
            A_forplot(:, :, i) = squeeze(A(i, :, :)).';
        end
        
        [A_tdomain_transpose] = istft(A_forplot, fs, Window=rectwin(NWIN_H), OverlapLength=NWIN_H-hopsize_H, FFTLength=NFFT_H, ConjugateSymmetric=true, FrequencyRange='onesided');
        A_tdomain = A_tdomain_transpose.';
         
        figure(15)
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

figure(16)
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

figure(17);
mesh(NumOfFrame_vector, freqs_vector, S_save_dB)
colorbar
clim([S_caxis_underlimit, S_caxis_upperlimit])
view(2)
xlim([1 NumOfFrame])
title('S\_save')
xlabel('frame')
ylabel('frequency(Hz)')

S_save_L_dB = mag2db(abs(S_save_L));

figure(18);
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
source_max  = max(abs(source(1, ((start_frame_A-1)*hopsize+NWIN):end)));
source_predict_max  = max(abs(source_predict(1, ((start_frame_A-1)*hopsize+NWIN):(Second-0.1)*fs)));
ratio_source_predict = source_max/source_predict_max;
source_predict = source_predict.*ratio_source_predict;

% 畫 source_predict time plot %
figure(19)
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
audiowrite('wav\source.wav', source(1, ((start_frame_A-1)*hopsize+NWIN):end), fs)

ratio_y = 1 / max(abs(y(look_mic, ((start_frame_A-1)*hopsize+NWIN):end))) ;
y_filemane_str = ['wav\y-', string(reverberation_time),'.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y(look_mic, ((start_frame_A-1)*hopsize+NWIN):end)*ratio_y, fs)

if strcmp(A_mode, 'H')
    source_predict_filemane_str = ['wav\source_predict-', string(reverberation_time), 'x', string(NWIN_H), 'x', string(NFFT_H), 'x', string(hopsize_H), 'x',...
        string(NWIN), 'x', string(NFFT), 'x',string(hopsize), '-A_mode=', A_mode, '.wav'];
elseif strcmp(A_mode, 'RssRsy')
    source_predict_filemane_str = ['wav\source_predict-', string(reverberation_time), 'x', string(NWIN_H), 'x', string(NFFT_H), 'x', string(hopsize_H), 'x',...
        string(NWIN), 'x', string(NFFT), 'x',string(hopsize), '-A_mode=', A_mode, '-RssRsy_ini_mode=', RssRsy_ini_mode, '.wav'];
end
source_predict_filemane = join(source_predict_filemane_str, '');
audiowrite(source_predict_filemane, source_predict(1, ((start_frame_A-1)*hopsize+NWIN):end), fs)
fprintf(source_predict_filemane)
fprintf('\n')

toc