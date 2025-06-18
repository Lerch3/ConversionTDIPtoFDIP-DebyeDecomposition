% Dimas Putra Rachmawan
% Modul
% Jurusan Teknik Geofisika
% Fakultas Teknik Sipil Perencanaan dan Kebumian
% Institut Teknologi Sepuluh Nopember
% Surabaya 2025
% Creating Reference Variances

clc; clear; close all;

%% --- KONFIGURASI AWAL ---
Nerrors = 1000;
relaxation_time = 0.5; % s
R0 = 1;

windows = logspace(-1, 0, 21); % 0.1 s to 1 s
t_ = sqrt(windows(1:end-1) .* windows(2:end));

f = @(t) 0.1 * exp(-t / relaxation_time); % Debye response

% Integrasi respons ideal
decay_ideal = zeros(size(t_));
for k = 1:length(t_)
    decay_ideal(k) = integral(f, windows(k), windows(k+1)) / (windows(k+1) - windows(k));
end

%% --- MODEL ERROR ---
error_model = @(x, a, b) x * a + b;
a = 1e-2;
b = 1e-6;
decay_var = (error_model(decay_ideal, a, b)).^2;

%% --- TAU GRID & FREKUENSI ---
tau_ = logspace(log10(min(t_)) - 1.5, log10(max(t_)) + 1.5, 100);
w = 2 * pi * 1;

%% --- INISIALISASI VARIABEL ---
estimates = zeros(Nerrors, 2); % Kolom 1: |Z|, Kolom 2: fase
errors = zeros(Nerrors, 1);
rmse_all = zeros(Nerrors, 1);
decays = zeros(Nerrors, length(t_));
logrtds = zeros(Nerrors, length(tau_));

%% --- SIMULASI ---
for i = 1:Nerrors
    rng(i); % reproducible noise
    decay_noisy = decay_ideal + randn(size(decay_ideal)) .* sqrt(decay_var);
    decays(i,:) = decay_noisy;

    % Inversi Debye (gunakan NNLS sederhana)
    A = exp(-t_' ./ tau_);
    m = lsqnonneg(A, decay_noisy(:));

    % Hitung Z kompleks
    Z_complex = R0 * (1 - sum(m(:)' ./ (1 + 1i * w * tau_)));
    estimates(i, 1) = abs(Z_complex);             % Magnitudo
    estimates(i, 2) = angle(Z_complex) * 1000;     % Fase dalam mrad

    % Estimasi RMSE
    rmse_all(i) = sqrt(mean((decay_noisy(:) - A * m).^2));

    % Estimasi error magnitude-phase dengan propagasi kovarian kasar
    Jz = [real(Z_complex), imag(Z_complex)] ./ (real(Z_complex)^2 + imag(Z_complex)^2);
    cov_Z = eye(2) * mean(decay_var); % kovarian kasar
    errors(i) = sqrt(Jz * cov_Z * Jz');

    logrtds(i,:) = log(m(:)' + 1e-10); % log tau domain
end

%% --- KOVARIAN ---
cov_M = cov(logrtds);
cov_E = cov_Z;
cov_ref = cov([log(estimates(:,1)), estimates(:,2)]);

%% --- VISUALISASI ---
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1)
plot(t_, decays(1:100:end,:)');
set(gca, 'XScale', 'log');
xlabel('t [s]'); ylabel('R_0\eta(t)');
title('(a) Noise realizations of the input transient'); grid on;

subplot(1,3,2)
histogram(logrtds(:,40), 'Normalization','pdf');
xlabel('m_{40}'); ylabel('Counts');
title('(b) Scatter of the sample model parameter m_{40}'); grid on;

subplot(1,3,3)
histogram(log(estimates(:,1)), 'FaceColor','g', 'Normalization','pdf'); hold on;
histogram(estimates(:,2), 'FaceColor','r', 'Normalization','pdf');
xlabel('log|Z| / \phi [mrad]'); ylabel('Counts');
title('(c) Scatter of the obtained FD estimates'); grid on;
legend('log|Z|', '\phi [mrad]');

%% --- KOVARIANSI ---
figure('Position', [100, 100, 1000, 300]);

subplot(1,3,1)
imagesc(cov_M); colorbar; title('(a) Covariance matrix C_M');
xlabel('m_j'); ylabel('m_i'); axis image;

subplot(1,3,2)
imagesc(cov_ref); colorbar; title('(b) Reference covariance matrix');
xlabel('param j'); ylabel('param i'); axis image;

subplot(1,3,3)
imagesc(cov_E); colorbar; title('(c) Covariance matrix C_E');
xlabel('Re/Im'); ylabel('Re/Im'); axis image;

fprintf('Rata-rata RMSE: %.4f\n', mean(rmse_all));
