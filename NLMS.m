function [e, v1_hat, W] = NLMS(d,v2,P,mu_hat,epsilon)
% Implementation of Normalized Least Mean Square adaptive filtering algorithm
% Usage: [e, v1_hat, W] = NLMS(d,v2,P,mu_hat,epsilon)

% Inputs:
% d: noisy ECG signal vector, measured by the primary sensor, shape:(1,number of samples)
% v2: pure noise vector, measured by the reference sensor, shape:(1, number of samples)
% P: the order of the filter
% mu_hat: the step size before division of NLMS algorithm
% epsilon: the small positive number added to the denominator

% Outputs:
% e: the estimated clean ECG signal vector, shape:(1,number of samples)
% v1_hat: the estimated noise vector, direct output of the filter, shape:(1,number of samples)
% W: the filter coefficients at each time step, shape:(P, number of samples)

% Initialization
num = length(d);
e = zeros(1, num);
w = zeros(P, 1);
W = zeros(P, num);
v1_hat = zeros(1, num);
v2_padding = [zeros(P-1,1);transpose(v2)]; % pad v2 for extraction of vector v2(n)

for n = 1:num
    v2_i = flip(v2_padding(n:n+P-1)); % get vector v2(n);[v2(n);v2(n-1);...v2(n-P+1)]
    v1_hat(n) = transpose(v2_i) * w; % estimate the noise
    e(n) = d(n) - v1_hat(n); % calculate the error(also the clean ECG signal)
    mu = mu_hat / (epsilon + norm(v2_i)^2); % calculate the time-varying step size
    w = w + mu * conj(v2_i) * e(n); % update the filter coefficients
    W(:,n) = w;
end
end
