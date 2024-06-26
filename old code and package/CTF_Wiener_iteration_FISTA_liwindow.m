clc; clear;
close all;

% 加入資料夾 %
addpath('wpe_v1.33')
addpath('cvx')
addpath('FISTA')

%% RIR parameter %%
SorNum = 1;                                              % source number
MicNum = 30;                                             % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)
Ts = 1/fs;                                               % Sample period (s)

% ULA %
MicStart = [1, 1.5, 1];
spacing = 0.02;
MicPos = zeros(MicNum, 3);
for i = 1:MicNum
    MicPos(i, :) = [MicStart(1, 1)+(i-1)*spacing MicStart(1, 2) MicStart(1, 3)];
end

SorPos = [2, 2.6, 1];                                    % source position (m)
room_dim = [5, 6, 2.5];                                  % Room dimensions [x y z] (m)
reverberation_time = 0.4;                                % Reverberation time (s)
points_rir = 8192;                                       % Number of rir points (需比 reverberation time 還長)
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 0;                                           % Disable high-pass filter

%% generate ground-truth RIR (h) %%
% load RIR 的 .mat 檔 %
rir_filename_str = ['h_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

% 畫 ground-truth RIR time plot %
figure(1)
look_mic = 10;
plot(h(look_mic, :), 'r');
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('h')
xlabel('points')
ylabel('amplitude')
shg

%% compute ground-truth CTF (H) %%
NFFT = 1024;
hopsize = 256;

% windows %
osfac = round(NFFT/hopsize);
win = hamming(NFFT);
coe = zeros(hopsize, 1);
for i = 1:osfac
    coe = coe + win((i-1)*hopsize+1:i*hopsize).^2;
end
coe = repmat(coe, [osfac, 1]);
awin = win./sqrt(NFFT*coe);    % analysis window
swin = win./sqrt(coe/NFFT);    % synthesis window

% compute H %
frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);
H = zeros(frequency, L, MicNum);
for k = 1:frequency
    zeta = xcorr(awin,swin).*exp(1i*2*pi*(k-1)*(-NFFT+1:NFFT-1)'/NFFT);
    for i = 1:MicNum
        H_temp = conv(h(i, :), zeta.');
        H(k, :, i) = H_temp(:, hopsize:hopsize:end);
    end
end

L_vector = 1:1:L;
freqs_vector = linspace(0, fs/2, frequency);
H_dB = mag2db(abs(H));

% 畫 ground-truth RIR frequency plot %
figure(2);
mesh(L_vector, freqs_vector, H_dB(:, :, look_mic))
colorbar
view(2)
xlim([1 L])
title('H')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% reconstruct RIR from H (h_reconstruct) %%
h_reconstruct = reconstruct_RIR(points_rir, NFFT, hopsize, L, win, frequency, MicNum, H);    % reconstruct_RIR(points_rir, NFFT, hopsize, L, win, frequency, MicNum, CTF)
 
ratio_h_reconstruct = max(abs(h(look_mic, :)))/max(abs(h_reconstruct(look_mic, :)));
h_reconstruct_forplot = h_reconstruct(look_mic, hopsize*(osfac-1)+1:end);

figure(3)
plot(h(look_mic, :), 'r');
hold on
plot(h_reconstruct_forplot*ratio_h_reconstruct, 'b');
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('reconstruct h from H')
legend('h', 'h\_reconstruct')
xlabel('points')
ylabel('amplitude')
shg

%% 讀音檔 or 產生 white noise source (source) %%
Second = 23;
SorLen =  Second*fs;

% load source %
[source_transpose, fs] = audioread('245.wav', [1, SorLen]);    % speech source
source = source_transpose.';

% source = wgn(1, SorLen, 0);    % white noise source

%% compute source signal for frequency (S) %%
% source 轉頻域 %
S = my_STFT(source, NFFT, hopsize, win);
NumOfFrame = size(S, 2);
NumOfFrame_vector = 1:1:NumOfFrame;
freqs_vector = linspace(0, fs/2, frequency);
S_dB = mag2db(abs(S));

% 畫 source frequency plot %
% figure(4);
% mesh(NumOfFrame_vector, freqs_vector, S_dB)
% colorbar
% view(2)
% xlim([1 NumOfFrame])
% title('S')
% xlabel('frame')
% ylabel('frequency(Hz)')
% shg

%% RIR mix source 先在時域上 convolution 再做 stft (y and Y) %%
% convolution source and RIR %
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(h(i, :), source);
end

extra_delay_y = (ceil(NFFT/hopsize) - 1)*hopsize;    % put delay for equilization between time convolution and CTF 
y = zeros(MicNum, SorLen);
y(:, extra_delay_y+1:end) = as(:, 1:SorLen-extra_delay_y);

% y 轉頻域 %
Y = zeros(frequency, NumOfFrame, MicNum);
for i = 1:MicNum
    Y(:, :, i) = my_STFT(y(i, :), NFFT, hopsize, win);
end

% 利用 CTF model 還原 mic frequency stft (Y_CTF) %
% Y_CTF = zeros(frequency, NumOfFrame, MicNum);
% for FrameNo= L:NumOfFrame
%     for i = 1 : MicNum
%         Y_CTF(:, FrameNo, i) = sum(H(:, :, i).*flip(S(:, FrameNo-L+1:FrameNo), 2), 2);    % frequency x 1 = sum(frequency x L .* frequency x L, 2)
%     end
% end
% 
% % Y_CTF 轉回時域 (y_ctf) %
% y_CTF = zeros(MicNum, SorLen);
% for i = 1:MicNum
%     y_CTF(i, :) = my_ISTFT(Y_CTF(:, :, i), hopsize, win);
% end
% 
% % 畫 mic signal 比較 time plot %
% ratio_y = 1 / max(abs(y(look_mic, :))) ;
% ratio_y_CTF = 1 / max(abs(y_CTF(look_mic, :))) ;
% 
% figure(5);
% plot(y(look_mic, :)*ratio_y, 'r');
% hold on
% plot(y_CTF(look_mic, :)*ratio_y_CTF, 'b');
% hold off
% ylim([-1.1 1.1])
% title('y and y\_ctf')
% xlabel('points')
% ylabel('magnitude')
% legend('y', 'y\_ctf')
% shg
% 
% % 存兩種 .wav 檔 %
% y_filemane_str = ['wav\y-', string(reverberation_time), '.wav'];
% y_filemane = join(y_filemane_str, '');
% audiowrite(y_filemane, y(look_mic, :)*ratio_y, fs)
% 
% y_CTF_filemane_str = ['wav\y_CTF-', string(reverberation_time), '.wav'];
% y_CTF_filemane = join(y_CTF_filemane_str, '');
% audiowrite(y_CTF_filemane, y_CTF(look_mic, :)*ratio_y_CTF, fs)

%% WPE %%
% y_nodelay = as(:, 1:SorLen);
% y_wpe = wpe(y_nodelay.', 'wpe_parameter.m');
% y_wpe = y_wpe.';
% 
% y_wpe_filename_str = ['y_all_wpe-', string(reverberation_time), '.mat'];
% y_wpe_filename = join(y_wpe_filename_str, '');
% save(y_wpe_filename, 'y_wpe')

y_wpe_filename_str = ['y_wpe-', string(reverberation_time), '.mat'];
y_wpe_filename = join(y_wpe_filename_str, '');
load(y_wpe_filename);

% ratio_y = 1 / max(abs(y_nodelay(look_mic, :))) ;
% ratio_y_wpe = 1 / max(abs(y_wpe(look_mic, :))) ;
% 
% figure(5);
% plot(y_nodelay(look_mic, :), 'r');
% hold on
% plot(y_wpe(look_mic, :), 'b');
% hold off
% ylim([-1.1 1.1])
% title('y and y\_wpe')
% xlabel('points')
% ylabel('magnitude')
% legend('y', 'y\_wpe')
% shg

%% DAS beamformer %%
a = zeros(MicNum, SorNum, frequency);
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((SorPos - MicPos(i, :)).^2));
end

for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, :, n) = exp(-1j*omega/c*distance)./distance;
end

w = a/MicNum;

Y_wpe = zeros(frequency, NumOfFrame, MicNum);
for i = 1:MicNum
    Y_wpe(:, :, i) = my_STFT(y_wpe(i, :), NFFT, hopsize, win);
end

Y_DAS = zeros(frequency, NumOfFrame);
for FrameNo= 1:NumOfFrame
    for n = 1: frequency
         Y_DAS(n, FrameNo) = w(:, :, n)'*squeeze(Y_wpe(n, FrameNo, :));
    end  

end

%% initial Rss Rsy %%
start_ini_frame = ceil((12*fs - NFFT)/hopsize) + 1;    % wpe 第12秒開始穩定
ini_frame = start_ini_frame+ 420;    % 一個 frame 約等於0.016秒  看要娶幾個來 initialize

% 矯正 Y_DAS 的 PSD 使之與 S 一樣 %
PSD_S = sum(abs(S(:, start_ini_frame:ini_frame)).^2, "all");
PSD_Y_DAS = sum(abs(Y_DAS(:, start_ini_frame:ini_frame)).^2, "all");
PSD_ratio = PSD_S/PSD_Y_DAS;
Y_DAS = Y_DAS*sqrt(PSD_ratio);

% 初始化 Rss %
forfac_ini_Rss = 0.999;
Rss = zeros(L, L, frequency);
for FrameNo= start_ini_frame:ini_frame
    S_choose = Y_DAS(:, FrameNo+1-L:FrameNo).';
    for n = 1: frequency
         Rss(:,:,n) = forfac_ini_Rss*Rss(:,:,n) + (1-forfac_ini_Rss)*flip(S_choose(:,n)) * flip(S_choose(:,n))';
    end  

end

% 初始化 Rsy %
forfac_ini_Rsy = 0.999;
Rsy = zeros(L, MicNum, frequency);
for FrameNo= start_ini_frame:ini_frame
    S_choose = Y_DAS(:, FrameNo+1-L:FrameNo).';
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
figure(6);
mesh(L_vector, freqs_vector, squeeze(A_ini_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A\_ini')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% A_ini 轉回時域 (A_ini_tdomain) %%
A_ini_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_ini_forplot(:, :, i) = squeeze(A_ini(i, :, :)).';
end

A_ini_tdomain = reconstruct_RIR(points_rir, NFFT, hopsize, L, win, frequency, MicNum, A_ini_forplot);
A_ini_tdomain_forplot = A_ini_tdomain(look_mic, hopsize*(osfac-1)+1:end);

figure(7)
plot(h(look_mic, :), 'r');
hold on
plot(A_ini_tdomain_forplot, 'b');
hold off
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('A\_ini\_tdomain')
legend('h', 'A\_ini\_tdomain')
xlabel('points')
ylabel('amplitude')
shg

%% initial A RAA rAy %%
A = zeros(MicNum, L, frequency);
RAA = zeros(L, L, frequency);
rAy = zeros(L, 1, frequency);

%% iteration with LASSO (cvx or FISTA) %%
alpha = 0.999;
beta = 0.999;
gamma = 0.001;
delta = 0.001;
lambda = 1e-4;
dia_load_A = 10^(-7);
dia_load_S_predict = 10^(-5);
A_mode = 'predict';    % 'predict' or 'true'
S_for_RssRsy = zeros(L, 1, frequency);
for FrameNo = ini_frame-(osfac-1)-L+1 : ini_frame-(osfac-1)
    S_for_RssRsy(ini_frame-(osfac-1)+1-FrameNo, :, :) = Y_DAS(:, FrameNo);
end

S_predict = zeros(L, 1, frequency);
S_save_first = zeros(frequency, NumOfFrame);
S_save_second = zeros(frequency, NumOfFrame);
S_save_third = zeros(frequency, NumOfFrame);
process_first_frame = ini_frame+1;

for FrameNo = process_first_frame : NumOfFrame
    for n = 1:frequency
        Y_now = squeeze(Y(n, FrameNo, :));
        Y_before = squeeze(Y(n, FrameNo-osfac, :));
        S_ini = squeeze(Y(n, FrameNo-L+1+(osfac-1):FrameNo+(osfac-1), 1)).';    % 拿第一顆麥克風的信號當 FISTA 的 initial
        if strcmp(A_mode, 'predict')
            % 更新 Rss Rsy %
            if FrameNo ~= process_first_frame
                Rss(:,:,n) = alpha.*Rss(:,:,n) + (1 - alpha).*S_for_RssRsy(:,:,n)*S_for_RssRsy(:,:,n)';
                Rsy(:,:,n) = beta.*Rsy(:,:,n) + (1 - beta).*S_for_RssRsy(:,:,n)*Y_before';    
            end

        end

        if strcmp(A_mode, 'predict')
            % 使用估到的 S 來更新 A %
            A(:, : ,n) = Rsy(:, :, n)'*inv(Rss(:, :, n) + dia_load_A.*eye(L));
        elseif strcmp(A_mode, 'true')
            % 使用 ground-truth H 來當作 A %
            if FrameNo == process_first_frame
                for i = 1 : MicNum
                    A(i, :, :) = H(:, :, i).';
                end

            end

        end

        % lasso for source prediction %
%         S_predict(:, :, n) = lasso_cvx(A(:, :, n), Y_now, lambda);
%         S_predict(:, :, n) = FISTA_CTF(Y_now, A(:, :, n), S_ini, lambda);
        S_predict(:, :, n) = FISTA_Xian(Y_now, A(:, :, n), lambda);
    end

    % 更新 S_for_RssRsy %
    if strcmp(A_mode, 'predict')
        S_for_RssRsy(2:end, :, :) = S_for_RssRsy(1:end-1, :, :);
        S_for_RssRsy(1, :, :) =  S_predict(osfac, :, :);
    end
       
    % 儲存 S_predict %
    S_save_first(:, FrameNo) = S_predict(osfac-1, :, :);
    S_save_second(:, FrameNo-(osfac-1)) = S_predict(osfac, :, :);
    S_save_third(:, FrameNo) = S_predict(osfac+1, :, :);

    % 畫每一個 frame 的 A_tdomain %
    A_forplot = zeros(frequency, L, MicNum);
    for i = 1 : MicNum
        A_forplot(:, :, i) = squeeze(A(i, :, :)).';
    end
    
    A_tdomain = reconstruct_RIR(points_rir, NFFT, hopsize, L, win, frequency, MicNum, A_forplot);
    A_tdomain_forplot = A_tdomain(look_mic, hopsize*(osfac-1)+1:end);
    
    figure(8)
    plot(h(look_mic, :), 'r');
    hold on
    plot(A_tdomain_forplot, 'b');
    hold off
    ylim([h_yaxis_underlimit h_yaxis_upperlimit])
    title('A\_tdomain', num2str(FrameNo))
    legend('h', 'A\_tdomain')
    xlabel('points')
    ylabel('amplitude')

    % print 出先在處理哪一個 frame %
    fprintf('right now processing frame = %d\n', FrameNo)
end

error_A_ini_tdomain = sum((A_ini_tdomain_forplot - h(look_mic, 1:points_rir-hopsize*(osfac-1))).^2);
error_A_tdomain = sum((A_tdomain_forplot - h(look_mic, 1:points_rir-hopsize*(osfac-1))).^2);

%% 畫圖看最後 A 在頻域的樣子 (A_fin) %%
A_dB = mag2db(abs(A));
figure(9);
mesh(L_vector, freqs_vector, squeeze(A_dB(look_mic, :, :)).')
colorbar
view(2)
xlim([1 L])
title('A\_fin')
xlabel('frame')
ylabel('frequency(Hz)')
shg

%% plot spectrogram (S_save) %%
S_save_dB = mag2db(abs(S_save_second));
figure(10);
mesh(NumOfFrame_vector, freqs_vector, S_save_dB)
colorbar
view(2)
xlim([process_first_frame NumOfFrame])
title('S\_save')
xlabel('frame')
ylabel('frequency(Hz)')

figure(11);
mesh(NumOfFrame_vector, freqs_vector, S_dB)
colorbar
view(2)
xlim([process_first_frame NumOfFrame])
title('S')
xlabel('frame')
ylabel('frequency(Hz)')

%% predicted source 還原回時域 (source_predict) %%
% S_save 轉回時域 %
source_predict_transpose = my_ISTFT(S_save_second, hopsize, win);
source_predict = source_predict_transpose.';

% adjust source_predict 的最高點使之與 source 的一樣 %
source_max  = max(abs(source(1, (ini_frame*hopsize+NFFT):end)));
source_predict_max  = max(abs(source_predict(1, (ini_frame*hopsize+NFFT):(Second-0.1)*fs)));
ratio_source_predict = source_max/source_predict_max;
source_predict = source_predict.*ratio_source_predict;

% 畫 source_predict time plot %
figure(12)
plot(source(1, :), 'r');
hold on
plot(source_predict(1, :), 'b');
hold off
title('source\_predict')
xlabel('points')
ylabel('magnitude')
legend('source', 'source\_predict')
shg

%% save .wav 檔 %%
audiowrite('wav\source_section.wav', source(1, (ini_frame*hopsize+NFFT):end), fs)

ratio_y = 0.8 / max(abs(y(look_mic, (ini_frame*hopsize+NFFT):end))) ;
y_filemane_str = ['wav\y_section-', string(reverberation_time), '.wav'];
y_filemane = join(y_filemane_str, '');
audiowrite(y_filemane, y(look_mic, (ini_frame*hopsize+NFFT):end)*ratio_y, fs)

source_predict_filemane_str = ['wav\source_predict_section_Wiener_FISTA-', string(reverberation_time), 'x', string(NFFT), 'x',string(hopsize), ...
    '-A_mode=', A_mode, '.wav'];
source_predict_filemane = join(source_predict_filemane_str, '');
audiowrite(source_predict_filemane, source_predict(1, (ini_frame*hopsize+NFFT):end), fs)

fprintf('done\n')