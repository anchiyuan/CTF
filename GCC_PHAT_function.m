function delay_time = GCC_PHAT_function(x1, x2, Fs)

%% resemple mic signals %%
fs = 16000;
x1 = resample(x1, fs, Fs);
x2 = resample(x2, fs, Fs);

%% Window parameters %%
NWIN = 1024;
hopsize = NWIN/2;
NumOfFrame = 2*floor(length(x1)/NWIN)-1;
win = hann(NWIN+1);
win = win(1:end-1).';

%% FFT %%
NFFT = NWIN;
sum = zeros(1, NFFT/2+1);
for FrameNo = 1:NumOfFrame
    t_start = (FrameNo-1)*hopsize;
    tt = (t_start+1):(t_start+NWIN);
    
    x1_win = x1(tt).*win;
    x1_zp = [x1_win, zeros(1, (NFFT-NWIN))];
    X1 = fft(x1_zp, NFFT);
    X1_half = X1(1:NFFT/2+1);
    
    x2_win = x2(tt).*win;
    x2_zp = [x2_win, zeros(1,(NFFT-NWIN))];
    X2 = fft(x2_zp, NFFT);
    X2_half = X2(1:NFFT/2+1);
    
    XX = conj(X1_half).*X2_half;

    sum = sum + XX;
    
end

G = sum/NumOfFrame;
G_out = G./abs(G);    % whitening
            
%% IFFT output %%
G_out(1, end) = abs(G_out(1, end));
G_total = [G_out, conj(fliplr((G_out(1, 2:end-1))))];
phi = ifft(G_total, NFFT);

%% cheak delay or advance %%
[~, s] = max(abs(phi));
ss = s;
if ss > (NWIN/2)
    s = NWIN + 1 - s;
end

%% Interpolation %%
tt = -NWIN/fs/2:1/(100*fs):NWIN/fs/2;
T = 1/(fs);
phi_t = zeros(length(tt), 1);
for m = 1:length(tt)
    phi_t(m) = 0;
    
    si = s-100;
    sf = s+100;
    if si < 1
        si = 1;
    end

    if sf > NWIN
        sf = NWIN;
    end
    
    for n = si:sf
        if tt(m) == n*T
            phi_t(m) = phi_t(m) + phi(n);
        else       
            phi_t(m) = phi_t(m) + phi(n)*sin(pi*(tt(m)-n*T)/T)/(pi*(tt(m)-n*T)/T);
        end

    end

end

[~,a] = max(abs(phi_t));
if ss > (NWIN/2)
    delay_time = -tt(a);
else
    delay_time = tt(a);
end
    