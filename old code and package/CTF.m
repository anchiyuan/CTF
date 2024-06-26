clc; clear; 
close all;

%% rir parameter %%
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
%     MicPos(mic_num, :) = [1 + radius_UCA*cosd((angle_step)*(mic_num - 1)) 3 + radius_UCA*sind((angle_step)*(mic_num - 1)) 1];    % Receiver position [x y z] (m)
% end

% ULA %
MicStart = [1, 1.5, 1];
spacing = 0.02;
MicPos = zeros(MicNum, 3);
for mic_num = 1 : MicNum
    MicPos(mic_num, :) = [MicStart(1,1) + mic_num*spacing MicStart(1,2) MicStart(1,3)];
end

SorPos = [2, 2.5, 1];
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
plot3( [0 room_dim(1,1) room_dim(1,1) 0 0 0 room_dim(1,1) room_dim(1,1) 0 0 room_dim(1,1) room_dim(1,1) 0 0 room_dim(1,1) room_dim(1,1)], ...
       [0 0 room_dim(1,2) room_dim(1,2) 0 0 0 room_dim(1,2) room_dim(1,2) room_dim(1,2) room_dim(1,2) room_dim(1,2) room_dim(1,2) 0 0 0], ...
       [0 0 0 0 0 room_dim(1,3) room_dim(1,3) room_dim(1,3) room_dim(1,3) 0 0 room_dim(1,3) room_dim(1,3) room_dim(1,3) room_dim(1,3) 0] , 'k')
hold on
plot3(MicPos(:,1), MicPos(:,2), MicPos(:,3), "*", 'MarkerSize', 3)
hold on
plot3(SorPos(:,1), SorPos(:,2), SorPos(:,3), ".", 'MarkerSize', 40)
hold off
xlabel('W = 5m')
ylabel('L = 6m')
zlabel('H = 2.5m','Rotation', 0, 'Position', [-0.7 -0.5 1])
shg

%% generate rir %%
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
look_mic = 15;
figure(2)
plot(h(look_mic,:), 'r');
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h')
xlabel('points')
ylabel('magnitude')
shg

%% 讀音檔 or 產生 white noise source %%
Second = 20;
SorLen =  Second*fs;

% load source %
[source, pathname] = read_audio(SorNum, Second, fs);     % audio source
% source = wgn(1, SorLen, 0);                              % white noise source

figure(3)
plot(source(1, :));
title('source')
xlabel('points')
ylabel('magnitude')
shg

%% compute ground-truth CTF (H) %%
NWIN = 64;
NFFT = 64;
hopsize = floor(NWIN/2);
L = floor((points_rir-NWIN)/hopsize) + 1;
L_pad = L+1;
NumOfFrame = floor((SorLen-NWIN)/hopsize) + 1;    % -1 是因為第一個 frame 長度是 hopsize 兩倍 ， floor 是因為最後一個 frame 不夠 hopsize 長度所以捨去
frequency = NWIN/2+1;
freqs_vec_half_fs = fs/2*linspace(0, 1, NWIN/2+1);

% RIR 轉頻域 %
h_transpose = h.';
h_pad = [zeros(MicNum,hopsize), h];
h_pad_transpose = h_pad.';
win_h = rectwin(NWIN).*ones(NWIN, MicNum);    % rectangular window
% win_h = hamming(NWIN).*ones(NWIN, MicNum);    % hamming window
win_h_pad = hamming(NWIN).*ones(NWIN, MicNum);
H = zeros(MicNum, L, frequency);
H_pad = zeros(MicNum, L_pad, frequency);
for FrameNo = 1 : L
    t_start = (FrameNo-1)*hopsize;
    tt = (t_start+1):(t_start+NWIN);
    h_win(tt, :) = h_transpose(tt, :).*win_h;
    H_temp = fft(h_win(tt, :), NWIN, 1);
    H_temp = H_temp(1:frequency, :).';
    H(:, FrameNo, :) = H_temp;

end

% for FrameNo = 1 : L_pad
%     t_start = (FrameNo-1)*hopsize;
%     tt = (t_start+1):(t_start+NWIN);
%     h_pad_win(tt, :) = h_pad_transpose(tt, :).*win_h_pad;
%     H_pad_temp = fft(h_pad_win(tt, :), NWIN, 1);
%     H_pad_temp = H_pad_temp(1:frequency, :).';
%     H_pad(:, FrameNo, :) =  H_pad_temp;
% 
% end

% 畫 ground-truth RIR frequency plot %
figure(3)
plot(freqs_vec_half_fs, 20*log10(abs(squeeze(H(1, :, :)))))
H_yaxis_underlimit = min(20*log10(abs(squeeze(H(1, end, :)))))-10;
ylim([H_yaxis_underlimit 0])
legend
title('H')
xlabel('frequency(Hz)')
ylabel('dB')
shg

[x_axis,y_axis] = meshgrid(1:L, freqs_vec_half_fs);
figure(23);
mesh(x_axis, y_axis, 20*log10(abs(squeeze(H(1, :, :)).')))
colorbar
caxis([H_yaxis_underlimit 0])
view(2)
xlim([0 L])
ylim([0 fs/2])
title('H')
xlabel('frame')
ylabel('frequency(Hz)')
shg
 
% figure(33)
% plot(freqs_vec_half_fs, 20*log10(abs(squeeze(H_pad(1, :, :)))))
% ylim([H_yaxis_underlimit 0])
% legend
% title('H\_pad')
% xlabel('frequency(Hz)')
% ylabel('dB')
% shg

%% reconstruct RIR from H or H_pad %%
for FrameNo = 1 : L
    H_sym = [squeeze(H(:, FrameNo, :)).'; conj(flipud(squeeze(H(:, FrameNo, 2:end-1)).'))];
    H_tdomain = ifft(H_sym, NWIN, 1, 'symmetric');
    t_start_overlap = (FrameNo-1)*hopsize;
    tt_overlap = (t_start_overlap+1):(t_start_overlap+NWIN);
    if FrameNo == 1
        h_recon(:, tt_overlap) = H_tdomain.';
        Overlap_term = H_tdomain((hopsize+1):end, :);
    else
        temp = [Overlap_term ; zeros(hopsize, MicNum)] + H_tdomain;    % 前一個segment的hopsize+1點到最後一點，與當前segment的整段長度相加
        h_recon(:, tt_overlap) = temp.';
        Overlap_term = temp((hopsize+1):end, :);                       % 保留原訊號長度NWIN的一半(hopsize)到end的距離，準備下一次疊加
    end

end

% for FrameNo = 1 : L_pad
%     H_pad_sym = [squeeze(H_pad(:, FrameNo, :)).'; conj(flipud(squeeze(H_pad(:, FrameNo, 2:end-1)).'))];
%     H_pad_tdomain = ifft(H_pad_sym, NWIN, 1, 'symmetric');
%     t_start_overlap = (FrameNo-1)*hopsize;
%     tt_overlap = (t_start_overlap+1):(t_start_overlap+NWIN);
%     if FrameNo == 1
%         h_pad_recon(:, tt_overlap) = H_pad_tdomain.';
%         Overlap_term = H_pad_tdomain((hopsize+1):end, :);
%     else
%         temp = [Overlap_term ; zeros(hopsize, MicNum)] + H_pad_tdomain;    % 前一個segment的hopsize+1點到最後一點，與當前segment的整段長度相加
%         h_pad_recon(:, tt_overlap) = temp.';
%         Overlap_term = temp((hopsize+1):end, :);                           % 保留原訊號長度NWIN的一半(hopsize)到end的距離，準備下一次疊加
%     end
% 
% end

% h_pad_recon = h_pad_recon(:, (hopsize+1):end);

% 畫 reconstruct RIR time plot %
figure(4)
plot(h_recon(1,:));
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h\_recon')
xlabel('points')
ylabel('magnitude')
shg

% figure(24)
% plot(h_pad_recon(1,:));
% ylim([h_yaxis_underlimit h_yaxis_upperlimit])
% title('h\_pad\_recon')
% xlabel('points')
% ylabel('magnitude')
% shg

%% 調整 H 大小以符合 ground-truth RIR %%
h_recon_max = max(h_recon(1,:));
h_max = max(h(1, :));
ratio = h_max/h_recon_max;
H = H.*ratio;

for FrameNo = 1 : L
    H_sym = [squeeze(H(:, FrameNo, :)).'; conj(flipud(squeeze(H(:, FrameNo, 2:end-1)).'))];
    H_tdomain = ifft(H_sym, NWIN, 1, 'symmetric');
    t_start_overlap = (FrameNo-1)*hopsize;
    tt_overlap = (t_start_overlap+1):(t_start_overlap+NWIN);
    if FrameNo == 1
        h_recon(:, tt_overlap) = H_tdomain.';
        Overlap_term = H_tdomain((hopsize+1):end, :);
    else
        temp = [Overlap_term ; zeros(hopsize, MicNum)] + H_tdomain;    % 前一個segment的hopsize+1點到最後一點，與當前segment的整段長度相加
        h_recon(:, tt_overlap) = temp.';
        Overlap_term = temp((hopsize+1):end, :);                       % 保留原訊號長度NWIN的一半(hopsize)到end的距離，準備下一次疊加
    end

end

figure(34)
plot(h_recon(1,:));
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h\_recon\_adjust')
xlabel('points')
ylabel('magnitude')
shg

%% compute source signal for frequency (S) %%
% source 轉頻域 %
source_transpose = source.';
win_source = hamming(NWIN).* ones(NWIN, SorNum);
S = zeros(SorNum, frequency, NumOfFrame);
for FrameNo=1:NumOfFrame
    t_start = (FrameNo-1)*hopsize;
    tt = (t_start+1):(t_start+NWIN);
    s_win(tt,:) = source_transpose(tt,:).*win_source;
    S_temp = fft(s_win(tt,:),NWIN, 1);
    S_temp = S_temp(1:frequency,:).';
    S(:, :, FrameNo) = S_temp;
end

%% RIR mix source 先在時域上處理再做 fft (Y) %%
% convolution source and RIR %
for i = 1:MicNum
    as(i,:) = conv(h(i,:),source);
end

as = as(:,1:SorLen);

% 加上 white noise 當作 interferer %
% SNR = 40;
% for ii = 1:MicNum
%     y(ii, :) = awgn(as(ii, :),SNR,'measured');
% end

% 不加 white noise %
y = as;

% y 轉頻域 %
y_transpose = y.';
win_y = hamming(NWIN).* ones(NWIN, MicNum);
Y = zeros(MicNum, frequency, NumOfFrame);
for FrameNo = 1:NumOfFrame
    t_start = (FrameNo-1)*hopsize;
    tt = (t_start+1):(t_start+NWIN);
    y_win(tt,:) = y_transpose(tt,:).*win_y;
    Y_temp = fft(y_win(tt,:),NWIN, 1);
    Y_temp = Y_temp(1:frequency, :).';
    Y(:, :, FrameNo) = Y_temp;
end

%% initial Rss rsy %%
ini_frame = floor(NumOfFrame/10);

% 初始化 Rss %
Rss = zeros(L, L, frequency);
for FrameNo= L:ini_frame
    S_choose = squeeze(S(:, :, FrameNo+1-L:FrameNo)).';
    for n = 1: frequency
         Rss(:,:,n) = Rss(:,:,n) + flipud(S_choose(:,n)) * flipud(S_choose(:,n))';
    end   
end

Rss = Rss ./ (ini_frame-L+1);

% 初始化 Rsy %
Rsy = zeros(L, MicNum, frequency);
for FrameNo= L:ini_frame
    S_choose = squeeze(S(:, :, FrameNo+1-L:FrameNo)).';
    Y_choose = squeeze(Y(:, :, FrameNo));
    for n = 1: frequency
         Rsy(:,:,n) = Rsy(:,:,n) + flipud(S_choose(:,n)) * Y_choose(:,n)';  
    end   
end

Rsy = Rsy ./ (ini_frame-L+1);

% 畫圖看初始 A 在時域的樣子 %
A_ini = zeros(MicNum, L, frequency);
dia_load_ini = 10^(-4);
for n = 1:frequency
    A_ini(:,:,n) = Rsy(:,:,n)'*inv(Rss(:,:,n) + dia_load_ini.*eye(L));
end

for FrameNo = 1 : L
    A_ini_sym = [squeeze(A_ini(:, FrameNo, :)).'; conj(flipud(squeeze(A_ini(:, FrameNo, 2:end-1)).'))]; 
    A_ini_tdomain = ifft(A_ini_sym, NWIN, 1, 'symmetric');
    t_start_overlap = (FrameNo-1)*hopsize;
    tt_overlap = (t_start_overlap+1):(t_start_overlap+NWIN);
    if FrameNo == 1
        ctf_ini_tdomain(:, tt_overlap) = A_ini_tdomain.';
        Overlap_term = A_ini_tdomain((hopsize+1):end, :);
    else
        temp = [Overlap_term ; zeros(hopsize, MicNum)] + A_ini_tdomain;    % 前一個segment的hopsize+1點到最後一點，與當前segment的整段長度相加
        ctf_ini_tdomain(:, tt_overlap) = temp.';
        Overlap_term = temp((hopsize+1):end, :);                           % 保留原訊號長度NWIN的一半(hopsize)到end的距離，準備下一次疊加
    end

end

figure(5)
plot(ctf_ini_tdomain(1, :));
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('ctf\_ini\_tdomain')
xlabel('points')
ylabel('magnitude')
shg

figure(25)
plot(freqs_vec_half_fs, 20*log10(abs(squeeze(A_ini(1, :, :)))))
ylim([H_yaxis_underlimit 0])
legend
title('A\_ini')
xlabel('frequency(Hz)')
ylabel('dB')
shg


figure(35);
mesh(x_axis, y_axis, 20*log10(abs(squeeze(A_ini(1, :, :)).')))
colorbar
caxis([H_yaxis_underlimit 0])
view(2)
xlim([0 L])
ylim([0 fs/2])
title('A\_ini')
xlabel('frame')
ylabel('frequency(Hz)')

ctf_ini_tdomain_filename_str = ['ctf_ini_tdomain_', 'NWIN=', string(NWIN), '.mat'];
ctf_ini_tdomain_filemane = join(ctf_ini_tdomain_filename_str, '');
save(ctf_ini_tdomain_filemane, 'ctf_ini_tdomain')

%% ctf_ini_tdomain 圖平移看看是否與h相似 %%
[max_h,max_h_idx] = max(h(1,:));
[max_ctf,max_ctf_idx] = max(ctf_ini_tdomain(1, :));
ctf_mod = zeros(1,4096);
ctf_mod(1, max_h_idx:end-(max_ctf_idx-max_h_idx)) = ctf_ini_tdomain(1, max_ctf_idx:end);
figure(6)
plot(ctf_mod(1, :), 'c');
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('ctf\_mod')
xlabel('points')
ylabel('magnitude')
shg

%% initial A %%
A = zeros(MicNum, L, frequency);

%% initial RAA rAy %%
RAA = zeros(L, L, frequency);
rAy = zeros(L, 1, frequency);

%% recursive process %%
alpha = 0.99;
beta = 0.99;
gamma = 0.001;
delta = 0.001;
dia_load_A = 10^(-2);
dia_load_S_predict = 10^(-2);
save_mode = 'front';    % 'front' or 'back'
S_predict = zeros(L, 1, frequency);
S_save = zeros(1, frequency, NumOfFrame);
S_save_L = zeros(1, frequency, NumOfFrame);
process_first_frame = ini_frame+1;

for FrameNo = process_first_frame : NumOfFrame
    for n = 1:frequency
        Y_before = squeeze(Y(:, n, FrameNo-1));
        Y_now = squeeze(Y(:, n, FrameNo));
        if FrameNo ~= ini_frame+1
            Rss(:,:,n) = alpha.*Rss(:,:,n) + (1 - alpha).*S_predict(:,:,n)*S_predict(:,:,n)';
            Rsy(:,:,n) = beta.*Rsy(:,:,n) + (1 - beta).*S_predict(:,:,n)*Y_before';
        end
        
%         A(:,:,n) = Rsy(:,:,n)'*inv(Rss(:,:,n) + dia_load_A.*eye(L));
        A = H;
        RAA(:,:,n) = gamma.*RAA(:,:,n) + (1 - gamma).*A(:,:,n)'*A(:,:,n);
        rAy(:,:,n) = delta.*rAy(:,:,n) + (1 - delta).*A(:,:,n)'*Y_now;
        S_predict(:,:,n) = inv(RAA(:,:,n) + dia_load_S_predict.*eye(L))*rAy(:,:,n);
    end

    if strcmp(save_mode,'front')
        % 存最前面的 frame %
        S_save(:, :, FrameNo) = S_predict(1,:,:);
        S_save_L(:, :, FrameNo) = S_predict(L,:,:);
    else
        % 存最後面的 frame %
        S_save(:, :, FrameNo-(L-1)) = S_predict(L,:,:);
        if FrameNo == NumOfFrame
            for count = 1:L
                S_save(:, :, NumOfFrame+1-count) = squeeze(S_predict(count,:,:)).';
            end
    
        end     
    end
    
    % 一直畫 A 的變化 %
%     for FrameNo_A = 1 : L
%         A_sym = [squeeze(A(:,FrameNo_A,:)).'; conj(flipud(squeeze(A(:, FrameNo_A, 2:end-1)).'))]; 
%         A_tdomain = ifft(A_sym, NWIN, 1, 'symmetric');
%         t_start_overlap = (FrameNo_A-1)*hopsize;
%         tt_overlap = (t_start_overlap+1):(t_start_overlap+NWIN);
%         if FrameNo_A == 1
%             ctf_tdomain(:, tt_overlap) = A_tdomain.';
%             Overlap_term = A_tdomain((hopsize+1):end, :);
%         else
%             temp = [Overlap_term ; zeros(hopsize, MicNum)] + A_tdomain;    % 前一個segment的hopsize+1點到最後一點，與當前segment的整段長度相加
%             ctf_tdomain(:, tt_overlap) = temp.';
%             Overlap_term = temp((hopsize+1):end, :);                           % 保留原訊號長度NWIN的一半(hopsize)到end的距離，準備下一次疊加
%         end
%     
%     end
%     
%     figure(7)
%     plot(ctf_tdomain(1,:));
%     ylim([h_yaxis_underlimit h_yaxis_upperlimit])
%     title('ctf\_tdomain', num2str(FrameNo))
%     xlabel('points')
%     ylabel('magnitude')
%     shg


%     fprintf('right now processing frame = %d\n', FrameNo)

end

%% 畫圖看最後 A 在時域的樣子 %%
for FrameNo = 1 : L
    A_sym = [squeeze(A(:,FrameNo,:)).'; conj(flipud(squeeze(A(:, FrameNo, 2:end-1)).'))]; 
    A_tdomain = ifft(A_sym, NWIN, 1, 'symmetric');
    t_start_overlap = (FrameNo-1)*hopsize;
    tt_overlap = (t_start_overlap+1):(t_start_overlap+NWIN);
    if FrameNo == 1
        ctf_tdomain(:, tt_overlap) = A_tdomain.';
        Overlap_term = A_tdomain((hopsize+1):end, :);
    else
        temp = [Overlap_term ; zeros(hopsize, MicNum)] + A_tdomain;    % 前一個segment的hopsize+1點到最後一點，與當前segment的整段長度相加
        ctf_tdomain(:, tt_overlap) = temp.';
        Overlap_term = temp((hopsize+1):end, :);                           % 保留原訊號長度NWIN的一半(hopsize)到end的距離，準備下一次疊加
    end

end

figure(7)
plot(ctf_tdomain(1,:), 'g');
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('ctf\_finish\_tdomain')
xlabel('points')
ylabel('magnitude')
shg

%% plot spectrogram %%
[x_axis,y_axis] = meshgrid(1:NumOfFrame,freqs_vec_half_fs);
figure(11);
mesh(x_axis,y_axis,20*log10(abs(squeeze(S))))
colorbar
caxis([-90, 40])
view(2)
xlim([0 NumOfFrame])
ylim([0 fs/2])
title('true source spectrogram(dB)')
xlabel('frame')
ylabel('frequency(Hz)')

% figure(12);
% mesh(x_axis,y_axis,20*log10(abs(squeeze(Y(1,:,:).*10))))
% colorbar
% caxis([-90, 40])
% view(2)
% xlim([0 NumOfFrame])
% ylim([0 fs/2])
% title('microphone received(dB)')
% xlabel('frame')
% ylabel('frequency(Hz)')

figure(13);
mesh(x_axis,y_axis,20*log10(abs(squeeze(S_save))))
colorbar
caxis([-90, 40])
view(2)
xlim([0 NumOfFrame])
ylim([0 fs/2])
title('predicted source spectrogram(dB)')
xlabel('frame')
ylabel('frequency(Hz)')

% figure(14);
% mesh(x_axis,y_axis,20*log10(abs(squeeze(S_save_L))))
% colorbar
% caxis([-90, 40])
% view(2)
% xlim([0 NumOfFrame])
% ylim([0 fs/2])
% title('predicted source spectrogram L(dB)')
% xlabel('frame')
% ylabel('frequency(Hz)')

%% predicted source 還原回時域 %%
source_predict = zeros(1, SorLen);
recover_first_frame = (ini_frame+1)-(L-1);
for FrameNo = recover_first_frame : NumOfFrame
    S_sym = [squeeze(S_save(1, :, FrameNo)).';conj(flipud(squeeze(S_save(1, 2:end-1, FrameNo)))).']; 
    S_tdomain = ifft(S_sym, NWIN, 1, 'symmetric');
    S_tdomain = S_tdomain.*hamming(NWIN);    % ifft 完蓋 window
    t_start_overlap = (FrameNo-1)*hopsize;
    tt_overlap = (t_start_overlap+1):(t_start_overlap+NWIN);
    if FrameNo == recover_first_frame
        source_predict(1, tt_overlap) = S_tdomain.';
        Overlap_term = S_tdomain((hopsize+1):end, 1);
    else
        temp = [Overlap_term;zeros(hopsize, 1)] + S_tdomain;    % 前一個segment的hopsize+1點到最後一點，與當前segment的整段長度相加
        source_predict(1, tt_overlap) = temp.';
        Overlap_term = temp((hopsize+1):end, 1);                % 保留原訊號長度NWIN的一半(hopsize)到end的距離，準備下一次疊加
    end

end

max_source_predict  = max(abs(source_predict));
max_source  = max(abs(source));
ratio_source_predict = max_source/max_source_predict;
source_predict = source_predict.*ratio_source_predict;

figure(8)
plot(source(1,:), 'r');
hold on
plot(source_predict(1,:), 'b');
hold off
ylim([-1 1])
title('source\_predict')
xlabel('points')
ylabel('magnitude')
legend('true source', 'predicted source')
shg

%% save .wav %%
audiowrite('source.wav', source, fs)
audiowrite('y.wav', y(1,:).*10, fs)
source_predict_filemane_str = ['source_predict_NWIN=',string(NWIN),'_L=', string(L), '.wav'];
source_predict_filemane = join(source_predict_filemane_str, '');
audiowrite(source_predict_filemane, source_predict, fs)
fprintf(source_predict_filemane)
fprintf('\n')
