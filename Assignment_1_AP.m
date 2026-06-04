clear; close all; clc;

% Set randomizer for reproducibility
rng(2);



%% Generate data based on the requirements from the assignment
% Constants
M = 5;
N = 20;
theta = [-20; 30];
f = [0.1; 0.3];
SNR = 20;
d = 2;
Delta = 1/2;

% Data generator function
[X,A,S] = gendata(M,N,Delta,theta,f,SNR);

% Plot the singular values of the generated data, separating noise from
% signals
figure;
plot(svd(X),'-o','LineWidth',1.5);
grid on;
xlabel('Index');
ylabel('Singular value');
title('Singular values: baseline case');

% Compare singular values under changes
cases = {
    'Baseline',              M,   N,   [-20;30],  [0.1;0.3];
    'Double N',              M,   2*N, [-20;30],  [0.1;0.3];
    'Double M',              2*M, N,   [-20;30],  [0.1;0.3];
    'Close angles',          M,   N,   [-20;-18], [0.1;0.3];
    'Close frequencies',     M,   N,   [-20;30],  [0.1;0.105];
};

figure;
hold on;
for c = 1:size(cases,1)
    label = cases{c,1};
    Mc = cases{c,2};
    Nc = cases{c,3};
    thetac = cases{c,4};
    fc = cases{c,5};

    Xc = gendata(Mc,Nc,Delta,thetac,fc,SNR);
    sv = svd(Xc);
    semilogy(1:length(sv), sv, 'o-', 'DisplayName', label, 'LineWidth', 1.3);
end
grid on;
xlabel('Index');
ylabel('Singular value');
legend('Location','best');
title('Singular values under different scenarios');

%% Noiseless / high-SNR correctness checks
% constants
SNR_high = 300;
M = 5;
N = 20;
theta = [-20; 30];
f = [0.1; 0.3];
d = 2;

[X,A,S] = gendata(M,N,Delta,theta,f,SNR_high);

% Recover theta and f using esprit algorithm, also use the joint estimation
% algorithm to recover both
theta_hat = esprit(X,d);
f_hat = espritfreq(X,d);
[theta_joint,f_joint] = joint(X,d,5);

% Match to ensure the estimates correspond to the correct counterpart
theta_hat = match_estimates(theta_hat, theta);
f_hat = match_estimates(f_hat, f);



fprintf('\n=== High-SNR correctness check ===\n');
fprintf('True theta:       %s\n', mat2str(theta.',4));
fprintf('ESPRIT theta:     %s\n', mat2str(theta_hat.',4));
fprintf('Joint theta:      %s\n', mat2str(theta_joint.',4));

fprintf('True f:           %s\n', mat2str(f.',4));
fprintf('ESPRIT f:         %s\n', mat2str(f_hat.',4));
fprintf('Joint f:          %s\n', mat2str(f_joint.',4));

%% Performance comparison versus SNR
% constants
M = 3;
N = 20;
theta_true = [-20; 30];
f_true = [0.1; 0.12];
d = 2;
m_smooth = 5;

SNRvec = 0:4:20;
Nruns = 1000;

theta_esprit_all = zeros(d,Nruns,length(SNRvec));
f_esprit_all = zeros(d,Nruns,length(SNRvec));

theta_joint_all = zeros(d,Nruns,length(SNRvec));
f_joint_all = zeros(d,Nruns,length(SNRvec));

for is = 1:length(SNRvec)
    SNR = SNRvec(is);

    for r = 1:Nruns
        X = gendata(M,N,Delta,theta_true,f_true,SNR);

        th_e = esprit(X,d);
        ff_e = espritfreq(X,d);
        [th_j,ff_j] = joint(X,d,m_smooth);

        th_e = match_estimates(th_e, theta_true);
        ff_e = match_estimates(ff_e, f_true);

        theta_esprit_all(:,r,is) = th_e;
        f_esprit_all(:,r,is) = ff_e;

        theta_joint_all(:,r,is) = th_j;
        f_joint_all(:,r,is) = ff_j;
    end

    fprintf('Finished SNR = %d dB\n', SNR);
end

% Means and standard deviations
theta_esprit_mean = squeeze(mean(theta_esprit_all,2));
theta_esprit_std  = squeeze(std(theta_esprit_all,0,2));

f_esprit_mean = squeeze(mean(f_esprit_all,2));
f_esprit_std  = squeeze(std(f_esprit_all,0,2));

theta_joint_mean = squeeze(mean(theta_joint_all,2));
theta_joint_std  = squeeze(std(theta_joint_all,0,2));

f_joint_mean = squeeze(mean(f_joint_all,2));
f_joint_std  = squeeze(std(f_joint_all,0,2));

% Plot angle mean
figure;
plot(SNRvec, theta_esprit_mean(1,:), 'o-', 'LineWidth',1.4); hold on;
plot(SNRvec, theta_esprit_mean(2,:), 'o-', 'LineWidth',1.4);
plot(SNRvec, theta_joint_mean(1,:), 'x--', 'LineWidth',1.4);
plot(SNRvec, theta_joint_mean(2,:), 'x--', 'LineWidth',1.4);
yline(theta_true(1),':');
yline(theta_true(2),':');
grid on;
xlabel('SNR [dB]');
ylabel('Mean estimated angle [deg]');
legend('ESPRIT theta 1','ESPRIT theta 2','Joint theta 1','Joint theta 2','Location','best');
title('Mean angle estimates versus SNR');

% Plot angle std
figure;
plot(SNRvec, theta_esprit_std(1,:), 'o-', 'LineWidth',1.4); hold on;
plot(SNRvec, theta_esprit_std(2,:), 'o-', 'LineWidth',1.4);
plot(SNRvec, theta_joint_std(1,:), 'x--', 'LineWidth',1.4);
plot(SNRvec, theta_joint_std(2,:), 'x--', 'LineWidth',1.4);
grid on;
xlabel('SNR [dB]');
ylabel('Std of angle estimates [deg]');
legend('ESPRIT theta 1','ESPRIT theta 2','Joint theta 1','Joint theta 2','Location','best');
title('Angle standard deviations versus SNR');

% Plot frequency mean
figure;
plot(SNRvec, f_esprit_mean(1,:), 'o-', 'LineWidth',1.4); hold on;
plot(SNRvec, f_esprit_mean(2,:), 'o-', 'LineWidth',1.4);
plot(SNRvec, f_joint_mean(1,:), 'x--', 'LineWidth',1.4);
plot(SNRvec, f_joint_mean(2,:), 'x--', 'LineWidth',1.4);
yline(f_true(1),':');
yline(f_true(2),':');
grid on;
xlabel('SNR [dB]');
ylabel('Mean estimated frequency');
legend('ESPRIT f 1','ESPRIT f 2','Joint f 1','Joint f 2','Location','best');
title('Mean frequency estimates versus SNR');

% Plot frequency std
figure;
plot(SNRvec, f_esprit_std(1,:), 'o-', 'LineWidth',1.4); hold on;
plot(SNRvec, f_esprit_std(2,:), 'o-', 'LineWidth',1.4);
plot(SNRvec, f_joint_std(1,:), 'x--', 'LineWidth',1.4);
plot(SNRvec, f_joint_std(2,:), 'x--', 'LineWidth',1.4);
grid on;
xlabel('SNR [dB]');
ylabel('Std of frequency estimates');
legend('ESPRIT f 1','ESPRIT f 2','Joint f 1','Joint f 2','Location','best');
title('Frequency standard deviations versus SNR');

%%  4. Zero-forcing beamformer correctness, no noise
% Constants
M = 3;
N = 20;
theta_true = [-20; 30];
f_true = [0.1; 0.12];
d = 2;

[X,A,S] = gendata(M,N,Delta,theta_true,f_true,SNR_high);

theta_hat = esprit(X,d);
f_hat = espritfreq(X,d);

theta_hat = match_estimates(theta_hat, theta_true);
f_hat = match_estimates(f_hat, f_true);

[W_theta,Ahat_theta] = zf_from_theta(M,Delta,theta_hat);
[W_freq,Ahat_freq,Shat] = zf_from_freq(X,f_hat);

Y_theta = W_theta' * X;
Y_freq = W_freq' * X;

fprintf('\n=== Zero-forcing check, high SNR ===\n');
fprintf('norm(W_theta^H X - S) / norm(S), order may differ: %.3e\n', ...
    norm(Y_theta - S,'fro')/norm(S,'fro'));
fprintf('norm(W_freq^H X - Shat) / norm(Shat): %.3e\n', ...
    norm(Y_freq - Shat,'fro')/norm(Shat,'fro'));

%%  5. Spatial responses at SNR = 10 dB
SNR = 10;
[X,A,S] = gendata(M,N,Delta,theta_true,f_true,SNR);

theta_hat = esprit(X,d);
f_hat = espritfreq(X,d)

theta_hat = match_estimates(theta_hat, theta_true);
f_hat = match_estimates(f_hat, f_true);

[W_theta,Ahat_theta] = zf_from_theta(M,Delta,theta_hat);
[W_freq,Ahat_freq,Shat] = zf_from_freq(X,f_hat);

theta_grid = -90:0.1:90;

figure;
for i = 1:d
    y = spatial_response(W_theta(:,i), M, Delta, theta_grid);
    plot(theta_grid,y,'LineWidth',1.4); hold on;
end
grid on;
xlabel('\theta [deg]');
ylabel('|w^H a(\theta)|');
title('Spatial responses: ZF beamformer from ESPRIT angle estimates');
xline(theta_true(1),':');
xline(theta_true(2),':');
legend('Beamformer 1','Beamformer 2','True DOA 1','True DOA 2','Location','best');

figure;
for i = 1:d
    y = spatial_response(W_freq(:,i), M, Delta, theta_grid);
    plot(theta_grid,y,'LineWidth',1.4); hold on;
end
grid on;
xlabel('\theta [deg]');
ylabel('|w^H a(\theta)|');
title('Spatial responses: ZF beamformer from ESPRIT frequency estimates');
xline(theta_true(1),':');
xline(theta_true(2),':');
legend('Beamformer 1','Beamformer 2','True DOA 1','True DOA 2','Location','best');


%% Bonus: minimax beamformer

Delta = 1/2;
theta_grid = -90:0.1:90;

[w1,t1,theta_grid,y1] = minimax_beamformer(M,Delta,theta_hat(1),theta_hat(2),theta_grid);
[w2,t2,theta_grid,y2] = minimax_beamformer(M,Delta,theta_hat(2),theta_hat(1),theta_grid);

W = [w1 w2];

figure;
plot(theta_grid,y1,'LineWidth',1.5); hold on;
plot(theta_grid,y2,'LineWidth',1.5);
grid on;
xlabel('\theta [deg]');
ylabel('|w^H a(\theta)|');
xline(theta_true(1),':');
xline(theta_true(2),':');
legend('Beamformer 1','Beamformer 2','True DOA 1','True DOA 2','Location','best');
title('Minimax zero-forcing beamformer spatial responses');


[W_freq,Ahat,Shat,theta_freq,t_freq,theta_grid,Y_freq] = minimax_from_freq(X,f_hat);

figure;
plot(theta_grid,Y_freq(1,:),'LineWidth',1.5); hold on;
plot(theta_grid,Y_freq(2,:),'LineWidth',1.5);
grid on;
xlabel('\theta [deg]');
ylabel('|w^H a(\theta)|');
xline(theta_true(1),':');
xline(theta_true(2),':');
legend('Beamformer 1','Beamformer 2','True DOA 1','True DOA 2','Location','best');
title('Frequency-based minimax zero-forcing beamformers');