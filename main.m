clear;
close all;
load("ECG_database.mat")

%% Original clean ECG signal
s1 = Org_Data1(2, :)/200;
s2 = Org_Data2(2, :)/200;
s3 = Org_Data3(2, :)/200;

%% Noise
FS = 500;
TIME = 10;
wn = 0.1 * wgn(1, 5000, 0, 1, 42);
% wn = wn;
bwn = Org_bw_noise(1, 1:5000)/200;
pli = 0.1*sin(2*pi*50*(1/FS:1/FS:TIME));  % 50Hz
man = Org_ma_noise(1, 1:5000)/200;
emn = Org_em_noise(1, 1:5000)/200;

%% Noisy ECG Signal
d1 = s1 + wn;
d2 = s1 + bwn;
d3 = s1 + pli;
d4 = s1 + man;
d5 = s1 + emn;

%% Apply LMS filter to estimate the clean ECG signal
P = 8;
mu_wn = 2/max(eig(wn'*wn));
mu_bwn = 2/max(eig(bwn'*bwn));
mu_pli = 2/max(eig(pli'*pli));
mu_man = 2/max(eig(man'*man));
mu_emn = 2/max(eig(emn'*emn));

[e1_wn, v_hat1_wn, W1_wn] = LMS(d1,wn,P,mu_wn);
[e1_bwn, v_hat1_bwn, W1_bwn] = LMS(d2,bwn,P,mu_bwn);
[e1_pli, v_hat1_pli, W1_pli] = LMS(d3,pli,P,mu_pli);
[e1_man, v_hat1_man, W1_man] = LMS(d4,man,P,mu_man);
[e1_emn, v_hat1_emn, W1_emn] = LMS(d5,emn,P,mu_emn);

%% Apply NLMS filter to estimate the clean ECG signal
mu_hat = 0.02;
P = 2;
epsilon = 0.05;
[e2_wn, v_hat2_wn, W2_wn] = NLMS(d1,wn,P,mu_hat,epsilon);
[e2_bwn, v_hat2_bwn, W2_bwn] = NLMS(d2,bwn,P,mu_hat,epsilon);
[e2_pli, v_hat2_pli, W2_pli] = NLMS(d3,pli,P,mu_hat,epsilon);
[e2_man, v_hat2_man, W2_man] = NLMS(d4,man,P,mu_hat,epsilon);
[e2_emn, v_hat2_emn, W2_emn] = NLMS(d5,emn,P,mu_hat,epsilon);

%% Apply RLS filter to estimate the clean ECG signal
P = 2;
lambda = 1;
[e3_wn, v_hat3_wn, W3_wn] = RLS(d1,wn,P,lambda);
[e3_bwn, v_hat3_bwn, W3_bwn] = RLS(d2,bwn,P,lambda);
[e3_pli, v_hat3_pli, W3_pli] = RLS(d3,pli,P,lambda);
[e3_man, v_hat3_man, W3_man] = RLS(d4,man,P,lambda);
[e3_emn, v_hat3_emn, W3_emn] = RLS(d5,emn,P,lambda);

%% plot original clean ECG signal
plot_original_ecg(s1)

%% plot ECG signal corrupted by noise, after LMS, after NLMS, and after RLS
Plot_signal(d1,e1_wn,e2_wn,e3_wn,'White Noise')
Plot_signal(d2,e1_bwn,e2_bwn,e3_bwn,'Baseline Wander Noise')
Plot_signal(d3,e1_pli,e2_pli,e3_pli,'PLI')
Plot_signal(d4,e1_man,e2_man,e3_man,'Muscle Artifact')
Plot_signal(d5,e1_emn,e2_emn,e3_emn,'Electrode Motion Artifact')

%% plot the error magnitude vs sample index of different algorithms
plot_error_magnitude(s1,e1_wn,e2_wn,e3_wn,'White Noise')
plot_error_magnitude(s1,e1_bwn,e2_bwn,e3_bwn,'Baseline Wander Noise')
plot_error_magnitude(s1,e1_pli,e2_pli,e3_pli,'PLI')
plot_error_magnitude(s1,e1_man,e2_man,e3_man,'Muscle Artifact')
plot_error_magnitude(s1,e1_emn,e2_emn,e3_emn,'Electrode Motion Artifact')

%% SNR Improvement LMS Data1
snr_wn_LMS = snr_improvement(s1, d1, e1_wn);
snr_bwn_LMS = snr_improvement(s1, d2, e1_bwn);
snr_pli_LMS = snr_improvement(s1, d3, e1_pli);
snr_man_LMS = snr_improvement(s1, d4, e1_man);
snr_emn_LMS = snr_improvement(s1, d5, e1_emn);

%% SNR Improvement NLMS Data1
snr_wn_NLMS = snr_improvement(s1, d1, e2_wn);
snr_bwn_NLMS = snr_improvement(s1, d2, e2_bwn);
snr_pli_NLMS = snr_improvement(s1, d3, e2_pli);
snr_man_NLMS = snr_improvement(s1, d4, e2_man);
snr_emn_NLMS = snr_improvement(s1, d5, e2_emn);

%% SNR Improvement RLS Data1
snr_wn_RLS = snr_improvement(s1, d1, e3_wn);
snr_bwn_RLS = snr_improvement(s1, d2, e3_bwn);
snr_pli_RLS = snr_improvement(s1, d3, e3_pli);
snr_man_RLS = snr_improvement(s1, d4, e3_man);
snr_emn_RLS = snr_improvement(s1, d5, e3_emn);

%% MSE LMS Data1
mse_wn_LMS = mse(s1, e1_wn);
mse_bwn_LMS = mse(s1, e1_bwn);
mse_pli_LMS = mse(s1, e1_pli);
mse_man_LMS = mse(s1, e1_man);
mse_emn_LMS = mse(s1, e1_emn);

%% MSE NLMS Data1
mse_wn_NLMS = mse(s1, e2_wn);
mse_bwn_NLMS = mse(s1, e2_bwn);
mse_pli_NLMS = mse(s1, e2_pli);
mse_man_NLMS = mse(s1, e2_man);
mse_emn_NLMS = mse(s1, e2_emn);

%% MSE RLS Data1
mse_wn_RLS = mse(s1, e3_wn);
mse_bwn_RLS = mse(s1, e3_bwn);
mse_pli_RLS = mse(s1, e3_pli);
mse_man_RLS = mse(s1, e3_man);
mse_emn_RLS = mse(s1, e3_emn);
