function S = FISTA_CTF(Y, A, Sini, lambda)

% Y : MicNum x 1
% A : MicNum x L
% Sini : L x 1
% lambda : regularization factor

maxit = 1000;    % max iteration

tao = power_ite(A, Sini);    % poximal gradient learning rate (using power iteration algorithm)

Ypow = sum(sum(abs(Y).^2));
objval = zeros(maxit,1);

Sold = Sini;
Z = Sold;
t = 1;
for it = 1:maxit
       
    objval(it) = sum(abs(A*Sold - Y).^2) + lambda*sum(abs(Sold));
    if it>1 && objval(it-1)-objval(it)<Ypow*1e-5
        break
    end
    
    D = Z-(A'*(A*Z - Y))/tao;    
    S = SoftThresh_prox(D, lambda/tao);
    
    oldt = t;
    t = (1+sqrt(1+4*t^2))/2;
    
    Z = S + (oldt-1)/t * (S-Sold);
    Sold = S;   
end

objval = objval(1:it);
