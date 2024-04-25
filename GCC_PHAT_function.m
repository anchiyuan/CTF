function delay_time=GCC_PHAT_function(x1,x2,Fs)
%%
% clear all
% clc
% close all
% 
% fs=8000;
% % sample=200
% % tau=sample/fs
% tau=0.025
% sample=tau*fs
% 
% [x1 Fs]=audioread('D:\畢業光碟_鍾群\research\Sound Sources\female_16k_5s.wav');
% x1=resample(x1,fs,Fs);  % sample rate 從Fs變成fs

%-------------------------% delay
% x2=delayseq(x1,tau,fs);
% % x1=( awgn( x1,30 ) ).';
% % x2=( awgn( x2,30 ) ).';
% x1=x1.';
% x2=x2.';
% e=length(x1)
% figure(1)
% plot(1:e,x1,1:e,x2)  


% -------------------------% advamce
% x2=x1;
% x1=delayseq(x1,tau,fs);
% % x1=( awgn( x1,20 ) ).';
% % x2=( awgn( x2,20 ) ).';
% x1=x1.';
% x2=x2.';

% e=length(x1)
% figure(1)
% plot(1:e,x1,1:e,x2)  

% -------------------------%
% figure(1)
% plot(t,x1,t,x2)
%% delay estimate
fs=8000;
x1=resample(x1,fs,Fs);
x2=resample(x2,fs,Fs);
%% Windowing
NWIN=256*2;
hopsize=NWIN/2;  %50% overlap
NumOfFrame=2*floor(length(x1)/NWIN)-1;
win = hann(NWIN+1);  %hanning window
win = win(1:end-1).';

%% FFT
NFFT=NWIN;
df=fs/NFFT;
Freqs=linspace(0,fs/2,NFFT/2+1);%0:df:(NFFT/2-1)*df;
sum=zeros(1,NFFT/2+1);
% tic
for FrameNo=1:NumOfFrame
    %time segment
    t_start = (FrameNo-1)*hopsize;
    tt = (t_start+1):(t_start+NWIN);
    
    x1_win = x1(tt).*win;  %windowing
    x1_zp = [x1_win zeros(1,(NFFT-NWIN))];  % Zero-padding
    X1 = fft(x1_zp,NFFT);
    X1_half = X1(1:NFFT/2+1);
    
    x2_win = x2(tt).*win;  %windowing
    x2_zp = [x2_win zeros(1,(NFFT-NWIN))];  % Zero-padding
    X2 = fft(x2_zp,NFFT);
    X2_half = X2(1:NFFT/2+1);
    
    XX=conj(X1_half).*X2_half;

    sum=sum+XX;
    
end

G=sum/NumOfFrame;

G_out=G./abs(G);        % whitening
            

    % IFFT 
G_out(1,end)=abs(G_out(1,end));
G_total=[G_out,conj(fliplr((G_out(1,2:end-1))))];
phi=ifft(G_total,NFFT);

% figure(2)
% plot((1:NWIN),abs(phi))
%% cheak delay or advance
[~,s]=max(abs(phi));
ss=s;
if ss > (NWIN/2)
    s=NWIN+1-s;
end
%% Interpolation
tt=-NWIN/fs/2:1/(10*fs):NWIN/fs/2;
T=1/(fs);
for m=1:length(tt)
    phi_t(m)=0;
    
    si=s-20;
    sf=s+20;
    if si<1
        si=1;
    end        
    if sf>NWIN
        sf=NWIN;
    end
    
    for n=si:sf
%     for n=1:length(phi)
        if tt(m)==n*T;
            phi_t(m)=phi_t(m)+phi(n);
        else       
            phi_t(m)=phi_t(m)+phi(n)*sin(pi*(tt(m)-n*T)/T)/(pi*(tt(m)-n*T)/T);
        end
    end
end

% toc 

% figure(3)
% plot(tt,abs(phi_t))

 [~,a]=max(abs(phi_t));
 if ss > (NWIN/2)
    delay_time=-tt(a);
 else
    delay_time=tt(a);
 end
end





    
            
            
            
            
            
            
            
            
            
            
            
            
            
    
    