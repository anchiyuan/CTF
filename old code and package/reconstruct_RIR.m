function atf = reconstruct_RIR(points_rir, NFFT, hopsize, L, win, frequency, MicNum, CTF)

ATF_frame = (points_rir - NFFT)/hopsize + 1;
impulse = zeros(1, ((ATF_frame+L-1)-1)*hopsize+NFFT);
impulse(:, (L-1)*hopsize + 1) = 1;
IMPULSE = my_STFT(impulse, NFFT, hopsize, win);

ATF = zeros(frequency, ATF_frame, MicNum);
for FrameNo= L:ATF_frame+L-1
    for i = 1 : MicNum
        ATF(:, FrameNo-L+1, i) = sum(CTF(:, :, i).*flip(IMPULSE(:, FrameNo-L+1:FrameNo), 2), 2);
    end
end

atf = zeros(MicNum, points_rir);
for i = 1:MicNum
    atf(i, :) = my_ISTFT(ATF(:, :, i), hopsize, win);
end