%% Hankel function
z = 1 + 1i; % complex input
H = besselh(0, 2, z);
disp(H);

%% Draw object domain and the locate the source 
% STEP 2
kb = 1;
lambda = 2 * pi / kb; % wavefield of the background field

% Rectangle domain (D) corners
x1 = 0; y1 = 0; % upper left corner
x2 = lambda; y2 = lambda; % lower right corner

% Location of the source
sx = 1/2 * lambda;
sy = 10 * lambda;

% Figure
figure; hold on; axis equal;
xlabel('x'); ylabel('y');
title('Figure 1: Configuration Sketch');

% Draw rectangle D
rectangle('Position', [x1, y1, x2 - x1, y2 - y1], 'EdgeColor', 'b', 'LineWidth', 2);
text(x2, -lambda/4, 'Domain D', 'HorizontalAlignment', 'center');

% Mark origin
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
text(0, -lambda/5, '(0,0)', 'HorizontalAlignment', 'right');

% Mark right corner
plot(lambda, lambda, 'ko', 'MarkerFaceColor', 'k');
text(lambda+ lambda/2, lambda+lambda/4, '({\lambda},{\lambda})', 'HorizontalAlignment', 'right');

% Plot source location
plot(sx, sy, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(sx, sy + lambda/4, 'Source s', 'HorizontalAlignment', 'center');

% Optional: draw a line from source to domain
%line([sx, x2], [sy, y2], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

% Adjust axis
%xlim([-lambda, sx + lambda]);
ylim([-sy/10, sy + 2*lambda]);
set(gca, 'YDir', 'reverse');  % because y increases downward
grid on

%% STEP 5
h = lambda / 20; % Grid spacing
x = 0:h:lambda; % x coordinates (21 points)
y = 0:h:lambda; % y coordinates (21 points)
Nx = length(x);
Ny = length(y);

[X, Y] = meshgrid(x, y);   % Y increases downwards by default in images
%% STEP 6
% Compute |ro - ro_s| at each grid point
R = sqrt((X - sx).^2 + (Y - sy).^2);

% Compute the incident field (2D array) using Hankel function of the second kind, order 0
u_inc = -1j/4 * besselh(0, 2, kb * R);  % H_0^(2)(kb * R)

% Plot real, imag and absolute values
figure;

% Real part
subplot(1,3,1);
imagesc(x, y, real(u_inc));
axis equal tight;
title('Re(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Imaginary part
subplot(1,3,2);
imagesc(x, y, imag(u_inc));
axis equal tight;
title('Im(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Magnitude
subplot(1,3,3);
imagesc(x, y, abs(u_inc));
axis equal tight;
title('|u^{inc}|');
xlabel('x'); ylabel('y');
colorbar;
% Reverse y-axis to match coordinate system orientation (y increases downward)
set(gca,'YDir','normal');

%% STEP 8
% Spatially varying wavenumber field k(p)
% Define contrast domain (a circular (or a suqare) object in the center of the domain)
% xc = lambda / 2;           % center x
% yc = lambda / 2;           % center y

xc = lambda / 4;          
yc = lambda / 4;           
radius = lambda / 6;       % radius of the object
%side = lambda / 3;         % side length of the square object

% Define wavenumber field k(p)
k_object = 1.5 * kb;       % Higher wavenumber inside the object
k_field = kb * ones(Ny, Nx);  % Initialize as background

% Set k(p) > kb inside the circular object
object_mask = ((X - xc).^2 + (Y - yc).^2) <= radius^2;

%%Set k(p) > kb inside the square object
%object_mask = abs(X - xc) <= side/2 & abs(Y - yc) <= side/2;

k_field(object_mask) = k_object;

% Define non-negative contrast function chi(ro)
chi = (k_field / kb).^2 - 1;

%% STEP 8
% % Smoothly varying kb inside the object (e.g., Gaussian variation)
% % Define contrast domain (a circular (or a suqare) object in the center of the domain)
% xc = lambda / 2;           % center x
% yc = lambda / 2;           % center y
% 
% radius = lambda / 6;       % radius of the object
% 
% % Define wavenumber field k(p)
% k_field = kb * ones(Ny, Nx);  % Initialize as background
% 
% % Set k(p) > kb inside the circular object
% object_mask = ((X - xc).^2 + (Y - yc).^2) <= radius^2;
% 
% dist_from_center = sqrt((X - xc).^2 + (Y - yc).^2);
% kb_variation = 0.5 * exp(-((dist_from_center / radius).^2));  % smoothly varies from 0.5 to 0
% k_field(object_mask) = kb + kb_variation(object_mask);    % variable kb inside object
% 
% % Define non-negative contrast function chi(ro)
% chi = (k_field / kb).^2 - 1;

%% STEP 9
% Plot
figure;
imagesc(x, y, chi);
axis equal tight;
colorbar;
title('\chi(\rho): Non-Negative Contrast Function');
xlabel('x'); ylabel('y');

%% STEP 10
% Define the receiver domain D^{rec} as a line segment
figure; hold on; axis equal;
xlabel('x'); ylabel('y');
title('Figure 1: Configuration Sketch');

% Draw rectangle D
rectangle('Position', [x1, y1, x2 - x1, y2 - y1], 'EdgeColor', 'b', 'LineWidth', 2);
text(x2, -lambda/4, 'Domain D', 'HorizontalAlignment', 'center');

% Mark origin
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
text(0, -lambda/5, '(0,0)', 'HorizontalAlignment', 'right');

% Mark right corner
plot(lambda, lambda, 'ko', 'MarkerFaceColor', 'k');
text(lambda+ lambda/1.5, lambda+lambda/6, '({\lambda},{\lambda})', 'HorizontalAlignment', 'right');

% Plot source location
plot(sx, sy, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(sx, sy + lambda/4, 'Source s', 'HorizontalAlignment', 'center');
% Mark source location
plot(sx, sy, 'ko', 'MarkerFaceColor', 'k');
text(sx+lambda/5, sy-lambda/3, '(1/2{\lambda},10{\lambda})', 'HorizontalAlignment', 'left');

% Define the reciever domain
x_D_rec_start = -lambda; % x-coordinate of the start point of receiver domain
x_D_rec_end = 2*lambda; % x-coordinate of the end point of receiver domain
y_D_rec = 1.5*lambda; % y-coordinate of the receiver domain (does not change)
% Draw the D^{rec} in Figure 1
line([x_D_rec_start, x_D_rec_end], [y_D_rec, y_D_rec], 'Color', [0 0 0], 'LineWidth', 3);
%text(2.1*lambda, y_D_rec, 'Receiver Domain D^{rec}', 'HorizontalAlignment', 'left');
text(0.5*lambda, y_D_rec+lambda/2, 'Receiver Domain D^{rec}', 'HorizontalAlignment', 'center');
% Mark left corner
plot(-lambda, 1.5*lambda, 'ko', 'MarkerFaceColor', 'k');
text(-1.1*lambda, y_D_rec, '(-{\lambda},1.5{\lambda})', 'HorizontalAlignment', 'right');
% Mark right corner
plot(2*lambda, 1.5*lambda, 'ko', 'MarkerFaceColor', 'k');
text(2.1*lambda, y_D_rec, '(2{\lambda},1.5{\lambda})', 'HorizontalAlignment', 'left');

% Adjust axis
ylim([-sy/10, sy + 2*lambda]);
set(gca, 'YDir', 'reverse');  % because y increases downward
grid on
%% STEP 11
% Introduce grid on D^{rec}
%x_D_rec_grids = x_D_rec_start : 30*h : x_D_rec_end;
x_D_rec_grids = x_D_rec_start : 3*h : x_D_rec_end;
M = length(x_D_rec_grids);

figure; hold on; axis equal;
xlabel('x'); ylabel('y');
title('Figure 1: Configuration Sketch');

% Draw rectangle D
rectangle('Position', [x1, y1, x2 - x1, y2 - y1], 'EdgeColor', 'b', 'LineWidth', 2);
text(x2, -lambda/4, 'Domain D', 'HorizontalAlignment', 'center');

% Mark origin
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
text(0, -lambda/5, '(0,0)', 'HorizontalAlignment', 'right');

% Mark right corner
plot(lambda, lambda, 'ko', 'MarkerFaceColor', 'k');
text(lambda+ lambda/6, lambda+lambda/6, '({\lambda},{\lambda})', 'HorizontalAlignment', 'right');

% Define the reciever domain
x_D_rec_start = -lambda; % x-coordinate of the start point of receiver domain
x_D_rec_end = 2*lambda; % x-coordinate of the end point of receiver domain
y_D_rec = 1.5*lambda; % y-coordinate of the receiver domain (does not change)
% Draw the D^{rec} in Figure 1
scatter(x_D_rec_start:3*h: x_D_rec_end, y_D_rec);
%text(2.1*lambda, y_D_rec, 'Receiver Domain D^{rec}', 'HorizontalAlignment', 'left');
text(0.5*lambda, y_D_rec+lambda/4, 'Receiver Domain D^{rec}', 'HorizontalAlignment', 'center');
% Mark left corner
plot(-lambda, 1.5*lambda, 'ko', 'MarkerFaceColor', 'k');
text(-1.1*lambda, y_D_rec, '(-{\lambda},1.5{\lambda})', 'HorizontalAlignment', 'right');
% Mark right corner
plot(2*lambda, 1.5*lambda, 'ko', 'MarkerFaceColor', 'k');
text(2.1*lambda, y_D_rec, '(2{\lambda},1.5{\lambda})', 'HorizontalAlignment', 'left');

% Adjust axis
ylim([-sy/10, sy + 2*lambda]);
set(gca, 'YDir', 'reverse');  % because y increases downward
grid on
%% STEP 13
N = size(chi,1) * size(chi,2);
chi_vec = reshape(chi, [N,1]); % contrast function chi(ro) vector
%% STEP 14
% % Build the system matrix A (M x N)
% % Define voxel area and midpoint values of the grids (voxels)
% area = h^2; % area of a voxel via uniformly grid spacing
% x_midpoint = x/2;
% y_midpoint = y/2;
% 
% [X_midpoint, Y_midpoint] = meshgrid(x_midpoint, y_midpoint);  
% X_midpoint = reshape(X_midpoint, [1, N]);
% Y_midpoint = reshape(Y_midpoint, [1, N]);
% 
% Green_dist = zeros(M,N);
% for m = 1:M  
%     % Compute |ro_m - ro_n'| at each grid point, ro_m -> receiver
%     % locations, ro_n' -> midpoint locations of the voxels
%     row_Green = sqrt((x_D_rec_grids(m) - X_midpoint).^2 + (y_D_rec - Y_midpoint).^2); 
%     Green_dist(m,:) = reshape(row_Green, [1,N]);
% end
% 
% % Compute the Green function (2D array) using Hankel function of the second kind, order 0
% Green_func = -1j/4 * besselh(0, 2, kb * Green_dist);  
% 
% 
% % Compute the incident field (2D array)  distance
% u_inc_dist = sqrt((X_midpoint - sx).^2 + (Y_midpoint - sy).^2); 
% 
% % Compute the incident field (2D array) using Hankel function of the second kind, order 0
% u_inc_discritized = -1j/4 * besselh(0, 2, kb * u_inc_dist);  % H_0^(2)(kb * R)
% 
% A = area * Green_func .* u_inc_discritized;



% Build the system matrix A (M x N)
% Define voxel area and midpoint values of the grids (voxels)
area = h^2; % area of a voxel via uniformly grid spacing
x_fine = 0:(h/2):lambda;      % 41 points: spacing = lambda/40
y_fine = 0:(h/2):lambda;
x_midpoint = x_fine(1:2:end);      % take every second value → 21 points
y_midpoint = y_fine(1:2:end);

[X_midpoint, Y_midpoint] = meshgrid(x_midpoint, y_midpoint);  
Green_dist = zeros(M,N);
for m = 1:M  
    % Compute |ro_m - ro_n'| at each grid point, ro_m -> receiver
    % locations, ro_n' -> midpoint locations of the voxels
    row_Green = sqrt((x_D_rec_grids(m) - X_midpoint).^2 + (y_D_rec - Y_midpoint).^2); 
    Green_dist(m,:) = reshape(row_Green, [1,N]);
end

% Compute the Green function (2D array) using Hankel function of the second kind, order 0
Green_func = -1j/4 * besselh(0, 2, kb * Green_dist); 
% Compute G_s matrix
G_s = kb^2 *area * Green_func;
% Vectorize 2D incident field
u_inc_vec = reshape(u_inc, [N,1]);
%Compute A matrix
A = G_s * diag(u_inc_vec);

rank(A, 1e-3)
%% STEP 15

sing_vals = svd(A);

% Plot singular values
figure() 
plot(1:length(sing_vals), sing_vals, 'x', 'MarkerSize', 10, 'LineWidth', 4); 
xlabel('Index', 'FontSize', 14);
ylabel('Singular Value', 'FontSize', 14);
title('Singular Values of A', 'FontSize', 16);
%set(gca, 'FontSize', 12);
set(gca,'YScale', 'log');
grid on;

%% STEP 16
% Compute the scattered field
u_sc = A * chi_vec;

% STEP 17
% Find minimum norm solution to contrast function vector

%[U_r, S_r, V_r] = svd(A);
chi_vec_mn = pinv(A) * u_sc;

%% STEP 18
% Make an image of the minimum norm solution
chi_mn = reshape(real(chi_vec_mn), [Nx, Ny]);
chi_mn = max(chi_mn, 0);
% Plot
figure;
imagesc(x, y, chi_mn);
axis equal tight;
colorbar;
title('\chi_{mn}(\rho): Minimum Norm Contrast Function');
xlabel('x'); ylabel('y');

%% STEP 19
% Repeat for different M

% Introduce new grid on D^{rec}
h_new = 6*h ; %new uniform grid
x_D_rec_grids_new = x_D_rec_start : h_new : x_D_rec_end;
M_new = length(x_D_rec_grids_new);

% Build the system matrix A (M_new x N)
% Define voxel area and midpoint values of the grids (voxels)
area_new = h_new^2; % area of a voxel via uniformly grid spacing

Green_dist_new = zeros(M_new,N);
for m = 1:M_new  
    % Compute |ro_m - ro_n'| at each grid point, ro_m -> receiver
    % locations, ro_n' -> midpoint locations of the voxels
    row_Green = sqrt((x_D_rec_grids_new(m) - X_midpoint).^2 + (y_D_rec - Y_midpoint).^2); 
    Green_dist_new(m,:) = reshape(row_Green, [1,N]);
end

% Compute the Green function (2D array) using Hankel function of the second kind, order 0
Green_func_new = -1j/4 * besselh(0, 2, kb * Green_dist_new);  

%A_new = area_new * Green_func_new .* u_inc_discritized; %u_inc_discritized does not depend on number of measurements (M)

G_s_new = area_new * Green_func_new;

A_new = G_s_new * diag(u_inc_vec); % u_inc_vec remains the same (independt of number of measurements)

rank(A_new, 1e-3)

sing_vals_new = svd(A_new);

% Plot singular values
figure()
plot(1:length(sing_vals_new), sing_vals_new, 'x', 'MarkerSize', 10, 'LineWidth', 4); 
xlabel('Index', 'FontSize', 14);
ylabel('Singular Value', 'FontSize', 14);
title('Singular Values of A', 'FontSize', 16);
%set(gca, 'FontSize', 12);
set(gca,'YScale', 'log');
grid on;

% Compute the scattered field
u_sc_new = A_new * chi_vec;


% Find minimum norm solution to contrast function vector
chi_vec_mn_new = pinv(A_new) * u_sc_new;

chi_mn_new = reshape(real(chi_vec_mn_new), [Nx, Ny]);
chi_mn_new = max(chi_mn_new, 0);
% Plot
figure;
imagesc(x, y, chi_mn_new);
axis equal tight;
colorbar;
title('\chi_{mn}(\rho): Minimum Norm Contrast Function');
xlabel('x'); ylabel('y');


%% STEP 20
% Add noise
%noise = rand(size(u_sc));
SNR_dB = 30; % e.g., 20 dB SNR
signal_power = norm(u_sc)^2 / length(u_sc);
noise_power = signal_power / (10^(SNR_dB/10));
noise = sqrt(noise_power) * randn(size(u_sc));  % Gaussian noise
u_sc_noisy = u_sc + noise;

chi_vec_mn_noisy = pinv(A) * u_sc_noisy;

chi_mn_noisy = reshape(real(chi_vec_mn_noisy), [Nx, Ny]);
chi_mn_noisy = max(chi_mn_noisy, 0);
% Plot
figure;
imagesc(x, y, chi_mn_noisy);
axis equal tight;
colorbar;
title('\chi_{mn}(\rho): Minimum Norm Contrast Function From Noisy Incident Field');
xlabel('x'); ylabel('y');

% Repeat for different M
%noise_new = rand(size(u_sc_new));
signal_power_new = norm(u_sc_new)^2 / length(u_sc_new);
noise_power_new = signal_power_new / (10^(SNR_dB/10));
noise_new = sqrt(noise_power_new) * randn(size(u_sc_new));  % Gaussian noise
u_sc_new_noisy = u_sc_new + noise_new;

chi_vec_mn_new_noisy = pinv(A_new) * u_sc_new_noisy;

chi_mn_new_noisy = reshape(real(chi_vec_mn_new_noisy), [Nx, Ny]);
chi_mn_new_noisy = max(chi_mn_new_noisy, 0);
% Plot
figure;
imagesc(x, y, chi_mn_new_noisy);
axis equal tight;
colorbar;
title('\chi_{mn}(\rho): Minimum Norm Contrast Function From Noisy Incident Field');
xlabel('x'); ylabel('y');
