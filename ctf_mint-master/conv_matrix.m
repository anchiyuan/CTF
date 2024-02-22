function G = conv_matrix(g,ifl)
% g: nch x fl
% ifl: length of inverse filter
%
% G: fl+ifl-1 x (ifl x nch)

[nch,fl] = size(g);

G = zeros(fl+ifl-1,ifl*nch);
for i = 1:ifl
    for j = 1:nch
        G(:,(j-1)*ifl+i) = [zeros(i-1,1);g(j,:).';zeros(ifl-i,1)];        
    end
end