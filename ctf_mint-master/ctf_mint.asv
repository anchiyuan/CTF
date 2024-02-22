function y = ctf_mint(x,rir,delta,md,nfft,shift)
% CTF MINT
%
% x: time-domain microphone signal, size (No. of microphones x signal length )
% rir: room impulse response, size (No. of microphones x No. of sources x rir length)
% delta: energy regularization factor
% md: modeling delay
% nfft: STFT window (frame) length
% shift: STFT frame shift
%
% y: estimated time-domain source signal, size (No. of sources x signal length)
%
% Author: Xiaofei Li, INRIA Grenoble Rhone-Alpes
% Copyright: Perception Team, INRIA Grenoble Rhone-Alpes
% The algorithm is described in the paper:
% Xiaofei Li, Sharon Gannot, Laurent Girin and Radu Horaud. Multisource MINT Using 
% the Convolutive Transfer Function. ICASSP, 2018.
%


if nargin<5
    nfft = 1024;
    shift = 256;
end
if nargin<4
    md = 4;
end
if nargin<3
    delta = 1e-3;
end

[I,J,Lr] = size(rir);

% windows
osfac = round(nfft/shift);
win = hamming(nfft);
coe = zeros(shift,1);
for i = 1:osfac
    coe = coe + win((i-1)*shift+1:i*shift).^2;
end
coe = repmat(coe,[osfac,1]);
awin = win./sqrt(nfft*coe);                       % analysis window
swin = win./sqrt(coe/nfft);                       % synthesis window

% STFT of microphone signals
X = my_STFT(x(1,:),nfft,shift,win);
for i = 2:I
    X(:,:,i) = my_STFT(x(i,:),nfft,shift,win);
end
Lx = size(X,2);

%% CTF MINT
La = length(shift:shift:Lr+2*nfft-2); % length of CTF  xcorr : 2*nfft - 1  rir : Lr

Lh = round((La-1)/(I/J-1));         % length of inverse filters
Ld = La+Lh-1;                       % length of target function

% Applying CTF MINT for each frequency
K = nfft/2+1;
Y = zeros(K,Lx,J);
for k =  1:K
    
    % CTF
    ak = zeros(I,J,La);
    zeta = xcorr(awin,swin).*exp(1i*2*pi*(k-1)*(-nfft+1:nfft-1)'/nfft);
    for i = 1:I
        for j = 1:J
            aijk = conv(squeeze(rir(i,j,:)),zeta);
            ak(i,j,:) = aijk(shift:shift:end);
        end
    end     
     
    % Desired response
    zet = zeta(shift:shift:end);
    d = zeros(Ld,1);  
    d(md+1:md+length(zet)) = zet;   
    
    % Convolution matrix
    A = zeros(Ld*J,Lh*I);
    for j = 1:J
        A(Ld*(j-1)+1:Ld*j,:) = conv_matrix(squeeze(ak(:,j,:)),Lh);
    end
    
    % Estimate inverse filters and apply inverse filtering
    for j = 1:J
        g = [zeros(Ld*(j-1),1);d;zeros(Ld*(J-j),1)];
        
        Delta = delta*norm(squeeze(ak(:,j,:)),'fro')^2;
        hmint = (A'*A+Delta*eye(Lh*I))\(A'*g);        
        hmint = reshape(hmint,[Lh,I]);
        
        Ykj = zeros(Lx+Lh-1,1);
        for i = 1:I
            Ykj = Ykj+conv(squeeze(X(k,:,i)).',hmint(:,i));
        end
        
        Y(k,:,j) = Ykj(1:Lx);
    end
end

y = my_ISTFT(Y(:,:,1),shift,win);
for j = 2:J
    y(:,j) = my_ISTFT(Y(:,:,j),shift,win);
end
y = y';