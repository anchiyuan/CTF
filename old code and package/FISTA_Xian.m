function S = FISTA_Xian(Y, A, lambda)

iterations=150;
t_iter = 5;
S = zeros(size(A, 2), 1);
S_old = S;
Z = S_old;

for times = 1:iterations 
    G = A'*(A*Z - Y);
   
    % Adagrad
    epsilon = 0.001;    % 小步
    eta = 0.01;         % learning rate
    t_last = 0.005;     % 從t=0.005開始搜索
    fitness_t_last = norm(A*(Z - t_last.*G) - Y);
    r = 0;
    for t_time = 1:t_iter
        fitness_t_next = norm(A*(Z - (t_last + epsilon).*G) - Y);
        gradient_ = (fitness_t_next - fitness_t_last)/epsilon;
        r = r + gradient_.^2;
        delta_theta = -eta/(epsilon + sqrt(r))*gradient_;
        
        t_last = t_last + delta_theta;
        fitness_t_last = fitness_t_next;
    end
    t = t_last;
    D = Z - t.*G;
   
    % Proximal
    gradient_D = D ./ (eps + abs(D)) ;    % eps = 2.2204e-16
    res = abs(D)- t*lambda;
    for ss= 1:size(S, 1)
        if res(ss) >= 0
            res(ss) = res(ss);
        else
            res(ss) = 0;
        end

    end

    S = res.*gradient_D;
    Z = S + (times-2)/(times+1)*(S - S_old);
    S_old = S;    
end
