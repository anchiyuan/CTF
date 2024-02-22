function error = obj_func_Kalman(L, x, start_ini_frame, ini_frame, Y_DAS, n, Y_delay, look_mic)
weight = zeros(L, 1);
K = zeros(L, 1);    % Kalman gain %
P = x(1)*eye(L);    % error covariance matrix %
R = x(2);    % measurement noise covariance matrix %
for FrameNo = start_ini_frame:ini_frame
    % no time update only have measurement update %
    K = P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).')*inv(conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P*flip(Y_DAS(n, FrameNo-L+1:FrameNo).') + R);
    weight = weight + K*(conj(Y_delay(n, FrameNo, look_mic)) - conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*weight);
    P = P - K*conj(flip(Y_DAS(n, FrameNo-L+1:FrameNo)))*P;
end

error = 0;
for FrameNo= start_ini_frame:ini_frame
    error = error+  abs(Y_delay(n, FrameNo, look_mic) - weight'*flip(Y_DAS(n, FrameNo-L+1:FrameNo).'));
end
