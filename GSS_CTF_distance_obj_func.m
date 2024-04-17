function cosine_similarity = GSS_CTF_distance_obj_func(dis, ang, MicNum, SorNum, MicPos, fs, c, NFFT, hopsize, points_rir, Y_wpe, Y_delay, mode)

%% parameters setting %%
win = hamming(NFFT);
osfac = round(NFFT/hopsize);
frequency = NFFT/2 + 1;
L = length(hopsize:hopsize:points_rir+2*NFFT-2);
freqs_vector = linspace(0, fs/2, frequency);
NumOfFrame = size(Y_wpe, 2);

%% DAS beamformer (Y_DAS) %%
% 算 mic 與 source 之距離 %
source_pos = [dis*sind(ang), dis*cosd(ang), 0];
distance = zeros(MicNum, SorNum);
for i = 1 : MicNum
    distance(i, :) =  sqrt(sum((source_pos - MicPos(i, :)).^2));
end

% 算 freefield ATF %
a = zeros(MicNum, SorNum, frequency);
for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    a(:, :, n) = exp(-1j*omega/c*distance)./distance;
end

% 算 DAS weight %
w = a/MicNum;

% 算 Y_DAS %
Y_DAS = zeros(frequency, NumOfFrame);
for FrameNo= 1:NumOfFrame
    for n = 1: frequency
         Y_DAS(n, FrameNo) = w(:, :, n)'*squeeze(Y_wpe(n, FrameNo, :));
    end  

end

%% predict CTF with Kalman filter (A) %%
start_ini_frame = L;
ini_frame = NumOfFrame;

% Kalman filter %
A = zeros(MicNum, L, frequency);
parfor i = 1:MicNum
    for n =1:frequency
        weight = zeros(L, 1);
        P = 0.5*eye(L);    % error covariance matrix %
        K = zeros(L, 1);    % Kalman gain %
        R = 10^(-3);    % measurement noise covariance matrix %
        for FrameNo = start_ini_frame:ini_frame
            % no time update only have measurement update %
            K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
            weight = weight + K*(conj(Y_delay(n, FrameNo, i)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
            P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;
        end
    
        A(i, :, n) = weight';
    end

end

%% A 轉回時域且算 NRMSPM (A_tdomain NRMSPM) %%
A_forplot = zeros(frequency, L, MicNum);
for i = 1 : MicNum
    A_forplot(:, :, i) = squeeze(A(i, :, :)).';
end

A_tdomain = reconstruct_RIR_normalwindow(points_rir, NFFT, hopsize, L, win, fs, frequency, MicNum, A_forplot);    % dimension = MicNum x (points_rir+(osfac-1)*hopsize)
A_tdomain = A_tdomain(:, hopsize*(osfac-1)+1:end);

ATF = fft(A_tdomain, points_rir, 2);
ATF = ATF(:, 1:frequency);

%% cosine similarity %%
RTF_from_CTF = ATF./ATF(1, :);
RTF_from_CTF = reshape(RTF_from_CTF.', [MicNum*frequency 1]);
ATF_from_freefield = zeros(MicNum, frequency);

if strcmp(mode, 'plane wave')
    kappa = [sind(ang), cosd(ang), 0];
    for n = 1:frequency
        omega = 2*pi*freqs_vector(n);
        ATF_from_freefield(:, n) = exp(1j*omega/c*kappa*MicPos.').';
    end

elseif strcmp(mode, 'spherical wave')
    ATF_from_freefield = squeeze(a);

end

RTF_from_freefield = ATF_from_freefield./ATF_from_freefield(1, :);
RTF_from_freefield = reshape(RTF_from_freefield.', [MicNum*frequency 1]);

cosine_similarity = abs(RTF_from_CTF'*RTF_from_freefield/norm(RTF_from_CTF)/norm(RTF_from_freefield));
