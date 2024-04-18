function s = GSS_RTF_distance_obj_func(dis, ang, NFFT, hopsize, fs, y_nodelay, MicNum, MicPos)

%% parameters setting %%
win = hamming(NFFT);
frequency = NFFT/2 + 1;
freqs_vector = linspace(0, fs/2, frequency);

%% tfestimate for RTF %%
RTF_tfestimate =  tfestimate(y_nodelay(1, :), y_nodelay.', win, NFFT-hopsize, NFFT).';    % frequency x MicNum

%% freefield point source model for RTF %%
source_pos = [dis*sind(ang), dis*cosd(ang), 0];

r = zeros(MicNum, 1);
for i = 1 : MicNum
    r(i, :) =  sqrt(sum((source_pos - MicPos(i, :)).^2));
end

ATF_freefield = zeros(MicNum, frequency);
for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    ATF_freefield(:, n) = exp(-1j*omega/c*r)./r;
end

RTF_freefield = (ATF_freefield./ATF_freefield(1, :)).';    % frequency x MicNum


%% objective function %%
RTF_tfestimate = reshape(RTF_tfestimate, [MicNum*frequency 1]);
RTF_freefield = reshape(RTF_freefield, [MicNum*frequency 1]);

gamma = RTF_tfestimate'*RTF_freefield/norm(RTF_tfestimate)/norm(RTF_freefield);
s = 1/(1-real(gamma));

