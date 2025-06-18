% Dimas Putra Rachmawan
% Modul
% Jurusan Teknik Geofisika
% Fakultas Teknik Sipil Perencanaan dan Kebumian
% Institut Teknologi Sepuluh Nopember
% Surabaya 2025
% Conversion TDIP to FDIP

clc; clear; close all;

%% --- INPUT DATA DARI EXCEL ---
filename = 'Data.xlsx';
sheet = 'Sheet';
data = readtable(filename, 'Sheet', sheet);

% Ambil lokasi, kedalaman, dan ρ atau resistivitas
x = data.X_Lokasi;
z = data.Depth;
rho0 = data.Resistivitas;

% Ambil η(t) atau chargeabilitas
eta = table2array(data(:, 4:end));

% Waktu gates dalam milidetik, ubah ke detik
t_gate = [0.02, 0.04, 0.08, 0.16] / 1000;

%% --- SETUP TAU LOG-SPASI & FITTING (DEBYE MODEL) ---
N_tau = 20;
tau = logspace(log10(0.005), log10(1), N_tau);  % relaksasi time constants

% Buat matrix eksponensial: A(i,j) = exp(-t_i / tau_j)
A = exp(-t_gate' ./ tau);   % ukuran: [n_time x N_tau]

% Lakukan inversi non-linear (least square): eta ≈ A * m
m_all = zeros(length(x), N_tau);
eta_fit = zeros(length(x), length(t_gate));

for i = 1:length(x)
    m = lsqnonneg(A, eta(i,:)');   % solusi least-square non-negatif
    m_all(i,:) = m';
    eta_fit(i,:) = (A * m)';       % hasil fit
end

%% --- HITUNG RESPON Z(ω), ρ, ϕ, σ″ PADA f = 1 Hz ---
f = 1; w = 2 * pi * f;
rho_mag = zeros(length(x), 1);
sigma_real = 1 ./ rho0;
phi_mrad = zeros(length(x), 1);
sigma_imag = zeros(length(x), 1);

for i = 1:length(x)
    m = m_all(i,:) * 0.1;  % scaling agar efek IP kecil dan fisik

    % Impedansi kompleks
    Z_complex = rho0(i) * (1 - sum(m ./ (1 + 1i * w * tau)));
    rho_mag(i) = abs(Z_complex);

    % Fase IP dari rasio σ'' / σ'
    phi_mrad(i) = atan(w * imag(1/Z_complex) / sigma_real(i));

    % Komponen imaginer konduktivitas dari m & tau (Hase-like)
    sigma_imag(i) = w * sum((m .* tau) ./ (1 + (w * tau).^2));
end

%% --- HITUNG PFE ---
f_dc = 0.1;  % frekuensi DC (0.1 Hz)
f_ac = 1;    % frekuensi AC (1 Hz)
log_factor = log10(f_ac / f_dc);
PFE = ((rho0 - rho_mag) ./ rho_mag) / log_factor;

%% --- HITUNG MF ---
MF = (PFE ./ rho0) * 2000;

%% --- BUAT GRID INTERPOLASI UNTUK PLOT ---
nx = 200; nz = 100;
xq = linspace(min(x), max(x), nx);
zq = linspace(min(z), max(z), nz);
[Xq, Zq] = meshgrid(xq, zq);

% Interpolasi
RHO_grid = griddata(x, z, rho_mag, Xq, Zq, 'natural');
PHI_grid = griddata(x, z, phi_mrad, Xq, Zq, 'natural');
sigma_imag(sigma_imag <= 0) = 1e-10;
SIG_grid = griddata(x, z, log10(sigma_imag), Xq, Zq, 'natural');
PFE_grid = griddata(x, z, PFE, Xq, Zq, 'natural');
MF_grid = griddata(x, z, MF, Xq, Zq, 'natural');

%% --- VISUALISASI TOMOGRAFI ---
figure('Position', [100, 100, 1000, 800]);

% Resistivitas
subplot(5,1,1)
contourf(Xq, Zq, RHO_grid, 30, 'LineColor', 'none');
colormap(gca, 'jet'); colorbar;
title('\rho pada 1 Hz [\Omega\cdotm]');
ylabel('Depth (m)'); set(gca, 'YDir', 'reverse');

% Fase IP
subplot(5,1,2)
contourf(Xq, Zq, PHI_grid, 30, 'LineColor', 'none');
colormap(gca, flipud(jet)); colorbar;
title('\phi pada 1 Hz [mrad]');
ylabel('Depth (m)'); set(gca, 'YDir', 'reverse');

% Sigma Imag
subplot(5,1,3)
contourf(Xq, Zq, SIG_grid, 30, 'LineColor', 'none');
colormap(gca, 'jet'); colorbar;
caxis([-2.5 -1.5]);
title('\sigma" pada 1 Hz [log_{10}(S/m)]');
xlabel('Distance (m)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');

% PFE
subplot(5,1,4)
contourf(Xq, Zq, PFE_grid, 30, 'LineColor', 'none');
colormap(gca, 'jet'); colorbar;
title('PFE pada 1 Hz [%]');
xlabel('Distance (m)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');

% MF
subplot(5,1,5)
contourf(Xq, Zq, MF_grid, 30, 'LineColor', 'none');
colormap(gca, 'jet'); colorbar;
title('MF pada 1 Hz');
xlabel('Distance (m)'); ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');