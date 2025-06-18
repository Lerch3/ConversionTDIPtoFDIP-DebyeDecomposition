% Dimas Putra Rachmawan
% Modul
% Jurusan Teknik Geofisika
% Fakultas Teknik Sipil Perencanaan dan Kebumian
% Institut Teknologi Sepuluh Nopember
% Surabaya 2025
% Checking Error Propagation

clc; clear; close all;

%% --- KONFIGURASI ---
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
estimates = zeros(Nerrors, 2); % [log|Z|, phi]
errors_mag_phase = zeros(Nerrors, 2);
rmse_all = zeros(Nerrors, 1);
decays = zeros(Nerrors, length(t_));
logrtds = zeros(Nerrors, length(tau_));

%% --- SIMULASI ---
for i = 1:Nerrors
    rng(i);
    decay_noisy = decay_ideal + randn(size(decay_ideal)) .* sqrt(decay_var);
    decays(i,:) = decay_noisy;

    % Inversi Debye (NNLS)
    A = exp(-t_' ./ tau_);
    m = lsqnonneg(A, decay_noisy(:));

    % Hitung impedansi kompleks
    Z = R0 * (1 - sum(m(:)' ./ (1 + 1i * w * tau_)));
    Zmag = abs(Z);
    Zphase = angle(Z) * 1000;
    estimates(i, :) = [log(Zmag), Zphase];

    % RMSE
    rmse_all(i) = sqrt(mean((decay_noisy(:) - A * m).^2));

    % Kovarian fase & magnitude
    Jfd = [real(Z), imag(Z)] / (real(Z)^2 + imag(Z)^2);
    C_Z = eye(2) * mean(decay_var);
    err = sqrt(Jfd * C_Z * Jfd');
    errors_mag_phase(i, :) = [err, err];

    % log respons time domain
    logrtds(i,:) = log(m(:)' + 1e-10);
end

%% --- KOVARIANSI ---
cov_estimated_mag = cov(errors_mag_phase(:,1));
cov_estimated_phase = cov(errors_mag_phase(:,2));
cov_ref = cov(estimates);

%% --- VISUALISASI ---
figure('Position', [100, 100, 1000, 300]);

subplot(1,2,1)
histogram(errors_mag_phase(:,1), 'Normalization','pdf'); hold on;
xline(mean(errors_mag_phase(:,1)), 'k', 'LineWidth', 1.5);
xline(std(estimates(:,1)), 'r--', 'LineWidth', 1.5);
xlabel('std(ln|Z|)'); ylabel('Counts');
title('(a) Error of ln|Z|'); legend('Analytical', 'Mean Analytical', 'Reference'); grid on;

subplot(1,2,2)
histogram(errors_mag_phase(:,2), 'Normalization','pdf'); hold on;
xline(mean(errors_mag_phase(:,2)), 'k', 'LineWidth', 1.5);
xline(std(estimates(:,2)), 'r--', 'LineWidth', 1.5);
xlabel('std(\phi) [mrad]'); ylabel('Counts');
title('(b) Error of Phase'); legend('Analytical', 'Mean Analytical', 'Reference'); grid on;

fprintf('Rata-rata RMSE: %.4f\n', mean(rmse_all));
