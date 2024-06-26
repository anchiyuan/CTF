function S = lasso_cvx(A,Y, lambda)

[~, n] = size(A);
cvx_begin quiet
    variable S(n) complex
    minimize( norm(A*S - Y, 2) + lambda*norm(S, 1))
cvx_end
