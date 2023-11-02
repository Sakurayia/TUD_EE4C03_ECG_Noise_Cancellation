function [e, v1_hat, W] = RLS(d,v2,P,lambda)
% Implementation of Recursive Least Squares adaptive filtering algorithm
% Usage: [e, W] = RLS(d,v2,P,lambda)

% Inputs:
% d: noisy ECG signal vector, measured by the primary sensor, shape:(1,number of samples)
% v2: pure noise vector, measured by the reference sensor, shape:(1, number of samples)
% P: the order of the filter
% lambda: the forgetting factor of exponentially weighted RLS algorithm

% Outputs:
% e: the estimated clean ECG signal vector, shape:(1,number of samples)
% v1_hat: the estimated noise vector, direct output of the filter, shape:(1,number of samples)
% W: the filter coefficients at each time step, shape:(P, number of samples)

% Initialization
num = length(d);    % get the number of samples
e = zeros(1, num);
w = zeros(P, 1);
W = zeros(P, num);
v1_hat = zeros(1, num);
v2_padding = [zeros(P-1,1);transpose(v2)]; % pad v2 for extraction of vector v2(n)
delta = 0.05;
p_inv = eye(P)/delta; % the initialization of inverse autocorrelation matrix

% Update for each time step
for n = 1:num
    v2_n = flip(v2_padding(n:n+P-1)); % get vector v2(n);[v2(n);v2(n-1);...v2(n-P+1)]
    z = p_inv * conj(v2_n);
    g = z / (lambda + transpose(v2_n) * z); % calculate the gain vector
    v1_hat(n) = transpose(v2_n) * w;   % estimate the noise
    e(n) = d(n) - v1_hat(n);   % calculate the error(also the clean ECG signal)
    w = w + e(n) * g;    % update the filter coefficients
    W(:,n) = w;
    p_inv = (p_inv - g * ctranspose(z)) / lambda; % update matrix p_inv
end
end

