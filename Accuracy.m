% Dimas Putra Rachmawan
% Modul
% Jurusan Teknik Geofisika
% Fakultas Teknik Sipil Perencanaan dan Kebumian
% Institut Teknologi Sepuluh Nopember
% Surabaya 2025
% Accuracy of the TD to FD Conversion

clc; clear; close all;

%% --- SETUP: RELAXATION TIME DAN GATE WINDOWS
relaxation_times = logspace(-2, 1, 30);  % τ_s dari 0.01–10 s
omega = 2 * pi * 1;  % frekuensi 1 Hz

phases = zeros(length(relaxation_times), 1);
phases_real = zeros(length(relaxation_times), 1);
lam_all = zeros(length(relaxation_times), 1);
rmse_all = zeros(length(relaxation_times), 1);

error_model = @(x) x * 1e-2 + 1e-6;  % model error gauss

for itime = 1:length(relaxation_times)
    tau_s = relaxation_times(itime);

    % Time gates: logspace untuk windows integrasi
    windows = logspace(-1, 0, 21);
    t = sqrt(windows(1:end-1) .* windows(2:end));  % log-mean

    % Debye decay ideal
    gamma = 0.1;
    f = @(t) gamma * exp(-t / tau_s);  % Debye decay
    decay_ideal = zeros(1, length(t));

    for j = 1:length(t)
        l = windows(j);
        r = windows(j+1);
        decay_ideal(j) = integral(f, l, r) / (r - l);
    end

    % Tambah noise ke data decay
    rng(100);
    noise = randn(size(decay_ideal)) .* error_model(decay_ideal);
    decay = decay_ideal + noise;

    % Buat grid tau
    tau_grid = logspace(log10(0.005), log10(1), 30);
    A = zeros(length(t), length(tau_grid));
    for m = 1:length(t)
        for n = 1:length(tau_grid)
            A(m, n) = exp(-t(m) / tau_grid(n));
        end
    end

    % Bobot W dari error model
    W = diag(1 ./ (error_model(decay_ideal).^2));
    lambda = 1e2;
    L = eye(length(tau_grid));

    % Inversi Tikhonov
    AtW = A' * W;
    m = (AtW * A + lambda^2 * L' * L) \ (AtW * decay');

    % Hitung Z kompleks dari m hasil inversi
    sum_response = sum((m(:) ./ (1 + 1i * omega * tau_grid(:))).');  % skalar
    Z_complex = 1 - sum_response;

    % Hitung Z ideal (ground truth)
    Z_ideal = gamma / (1 + 1i * omega * tau_s);

    % Fase hasil inversi vs ideal
    phases(itime) = atan2(imag(Z_complex), real(Z_complex)) * 40;  % mrad
    phases_real(itime) = atan2(imag(Z_ideal), real(Z_ideal)) * 40;

    % Simpan lambda dan RMSE
    lam_all(itime) = lambda;
    rmse_all(itime) = sqrt(mean((A * m - decay').^2));
end

%% --- PLOT FIGURE (Mirip Gambar di Jurnal Hase)
figure('Position', [100 100 1000 400]);

% (a) Fase
subplot(1,2,1)
semilogx(relaxation_times, -phases_real, 'k.', 'DisplayName', 'Expected \phi');
hold on;
semilogx(relaxation_times, -phases, 'k-', 'DisplayName', 'Estimated \phi');
xline(windows(1), 'r:'); xline(windows(end), 'r:');
grid on;
xlabel('\tau_s [s]');
ylabel('-\phi [mrad]');
legend('Expected \phi', 'Estimated \phi');
title('(a) Phase accuracy');

% (b) Lambda dan RMSE
subplot(1,2,2)
yyaxis left
loglog(relaxation_times, lam_all, 'k-');
ylabel('\lambda_{final}');
yyaxis right
semilogx(relaxation_times, rmse_all, 'k--');
ylabel('\epsilon (RMSE)');
xlabel('\tau_s [s]');
xline(windows(1), 'r:'); xline(windows(end), 'r:');
legend('\lambda', '\epsilon');
grid on;
title('(b) Regularization and error');
