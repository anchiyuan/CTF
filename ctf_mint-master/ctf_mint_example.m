clc; clear; 
close all;

fs = 16000;

load 'rir.mat'
load 'source.mat'

[I,J,Lr] = size(rir);
[~,Ls] = size(source);

Lx = Ls+Lr-1;
x = zeros(I,Lx);
for i = 1:I
    for j = 1:J
        x(i,:) = x(i,:)+conv(source(j,:),squeeze(rir(i,j,:)));
    end
end

y = ctf_mint(x,rir,1e-5,6);

figure(1);
plot(source(2, :)', 'r')
hold on
plot(y(2, :)', 'b')
hold off
legend('source', 'x')

audiowrite('y.wav', y(2, :), fs)
audiowrite('source.wav', source(2, :), fs)

figure(2);
plot(squeeze(rir(3,1,:)))
title('RIR')