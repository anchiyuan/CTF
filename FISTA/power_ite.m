function tao = power_ite(A, Sini)

maxitl = 1000;

Snorm = sqrt(sum(abs(Sini).^2));
Sunit = Sini/Snorm;

for it=1:maxitl
    Sit = A'*A*Sunit;
    Snorm = sqrt(sum(abs(Sit).^2));
    
    if abs(sum(Sunit.*conj(Sit/Snorm)))>0.999
        break;
    else
        Sunit = Sit/Snorm;
    end
end
  
tao = 1.01*Snorm;
