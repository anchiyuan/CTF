clc; clear;
close all;

tic

%% RIR parameter %%
SorNum = 1;                                              % source number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MicNum = 9;                                             % number of microphone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)

mic_x = [ 100 ; 400 ; 400 ; 100 ; 100 ; 400 ; 400 ; 100 ; 250 ]./100;
mic_y = [ 100 ; 100 ; 400 ; 400 ; 100 ; 100 ; 400 ; 400 ; 250 ]./100;
mic_z = [ 50  ;  50 ;  50 ;  50 ; 200 ; 200 ; 200 ; 200 ; 200 ]./100;

MicPos = [mic_x, mic_y, mic_z,];


SorPos = [3.5, 1.7, 1];                                    % source position (m)
room_dim = [5, 6, 2.5];                                  % Room dimensions [x y z] (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reverberation_time = 0.2;                                % Reverberation time (s)
points_rir = 2048;                                       % Number of rir points (需比 reverberation time 還長)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 1;                                           % Disable high-pass filter

% 畫空間圖 %
figure(1);
plot3( [0 room_dim(1, 1) room_dim(1, 1) 0 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1)], ...
       [0 0 room_dim(1, 2) room_dim(1, 2) 0 0 0 room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) 0 0 0], ...
       [0 0 0 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0] , 'k')
hold on
plot3(MicPos(:, 1), MicPos(:, 2), MicPos(:, 3), 'r.', 'MarkerSize', 15)
hold on
plot3(SorPos(:, 1), SorPos(:, 2), SorPos(:, 3), '*', 'MarkerSize', 20)
hold off
xlabel('x\_axis')
ylabel('y\_axis')
zlabel('z\_axis')
title('空間圖')
shg

referencce_point = MicPos(1, :);
sorpos_groundtruth = SorPos - referencce_point;

%% generate ground-truth RIR (h) %%
% 產生 RIR 和存.mat 檔 %
h = rir_generator(c, fs, MicPos, SorPos, room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
rir_filename_str = ['h\h_for_localization_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
save(rir_filemane, 'h')

% load RIR 的 .mat 檔 %
rir_filename_str = ['h\h_for_localization_', string(reverberation_time), 'x', string(MicNum), 'x', string(points_rir), '.mat'];
rir_filemane = join(rir_filename_str, '');
load(rir_filemane)

look_mic = 1;
h_yaxis_upperlimit = max(h(look_mic, :)) + 0.01;
h_yaxis_underlimit = min(h(look_mic, :)) - 0.01;
% 畫 ground-truth RIR time plot %
figure(2)
plot(h(look_mic, :), 'r');
ylim([h_yaxis_underlimit h_yaxis_upperlimit])
title('ground-truth RIR')
xlabel('points')
ylabel('amplitude')
shg

[~, argmax_tf] = max(abs(h.'));

%% window parameter %%
NFFT = 1024;
hopsize = 256;
win = hamming(NFFT);
frequency = NFFT/2 + 1;
freqs_vector = linspace(0, fs/2, frequency);

%% 讀音檔 or 產生 white noise source (source) %%
Second = 23;
SorLen =  Second*fs;

% load speech source %
% [source_transpose, fs] = audioread('245.wav', [1, SorLen]);    % speech source
% source = source_transpose.';

% load white noise source %
source = wgn(1, SorLen, 0);                                    % white noise source

%% 產生麥克風訊號，在時域上 (y) %%
% convolution source and RIR %
as = zeros(MicNum, points_rir+SorLen-1);
for i = 1 : MicNum
    as(i, :) = conv(h(i, :), source);
end

y = as(:, 1:SorLen);

%% TDOA localization %%
% GCC-PHAT for TDOA estimation %
delay = zeros(MicNum-1, 1);
distance = zeros(MicNum-1, 1);
for i = 1:MicNum-1
    delay(i, :) = GCC_PHAT_function(y(1,:), y(i+1,:), fs);
    distance(i, :) = delay(i, :)*c;
end

MicPos = MicPos - referencce_point;
distance_groundtruth = zeros(MicNum-1, 1);
for i = 1:MicNum-1
    distance_groundtruth(i, :) = sqrt(sum((sorpos_groundtruth - MicPos(i+1, :)).^2)) - sqrt(sum((sorpos_groundtruth - MicPos(1, :)).^2));
end

distance = distance_groundtruth;

A = 2*[MicPos(2:end, :), distance];
b = MicPos(2:end, 1).^2 + MicPos(2:end, 2).^2 + MicPos(2:end, 3).^2 - distance.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = 4*sigma^2*(ones(MicNum-1, 1)*ones(MicNum-1, 1).' + eye(MicNum-1))*distance*distance.';
dia_load_psi = 1;
invpsi = inv(psi+dia_load_psi*eye(MicNum-1)); 
P = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 -1];

syms L
F = pinv( A.'*invpsi*A+L*P )*( A.'*invpsi*b );
c1 = F.'*P*F;
[num den]=numden(c1);  %把分子分母多項式分離出來
num2 = sym2poly(num);   %把多項式係數分離出來，變成向量式
RR=roots(num2);          %解所有的根(lamda)
real=RR( find(imag(RR)==0) );   %提出只有實數的根(lamda)
if length(real)==0 %|| rank(Psi)<(M-1)
     real=[real;0];
end

% 把解出的根帶入cos function ，當cos function最小時所對應的根(lamda)就是所求
if size(real,1)>1
    x = pinv( A.'*invpsi*A+real(1)*P )*( A.'*invpsi*b);
    cosfun1=(A*x-b).'*invpsi*(A*x-b);
    L2=real(1);
    for ii=2:length(real)
        x = pinv( A.'*invpsi*A+real(ii)*P )*( A.'*invpsi*b );
        cosfun2=(A*x-b).'*invpsi*(A*x-b);
        if cosfun2<cosfun1
            L2=real(ii);
            cosfun1=cosfun2;
        end
    end
else
    L2=real;
end
L2;

sorpos_estimation = ( pinv( A.'*invpsi*A+L2*P )*( A.'*invpsi*b ) );

toc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% threshold = 1e-2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iteration_time = 0;
% lambda_old = -1;
% lambda_new = 1;
% theta_old = pinv(A.'*invpsi*A+lambda_old*P)*A.'*invpsi*b;
% theta_new = pinv(A.'*invpsi*A+lambda_new*P)*A.'*invpsi*b;
% I_old = theta_old.'*P*theta_old;
% I_new = theta_new.'*P*theta_new;
% while 1
%     iteration_time = iteration_time + 1;
%     lambda_tmp = lambda_new;
%     I_tmp = I_new;
%     lambda_new = lambda_old - I_new*(lambda_new-lambda_old)/(I_new-I_old);
%     if isnan(lambda_new)
%         break;
%     end
%     theta_new = pinv(A.'*invpsi*A+lambda_new*P)*A.'*invpsi*b;
%     I_new = theta_new.'*P*theta_new;
%     if abs(lambda_new-lambda_tmp) < threshold
%         sorpos_estimation = theta_new(1:3, :).';
%         break;
%     end
% 
%     lambda_old = lambda_tmp;
%     I_old = I_tmp;
% end