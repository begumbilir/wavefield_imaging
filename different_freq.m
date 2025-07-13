%% STEP 21

%% Draw object domain and the locate the source 
% STEP 2
kb_0 = 1;
lambda = 2 * pi / kb_0; % wavefield of the background field

% Rectangle domain (D) corners
x1 = 0; y1 = 0; % upper left corner
x2 = lambda; y2 = lambda; % lower right corner

% Location of the source
sx = 1/2 * lambda;
sy = 10 * lambda;

% % Figure
% figure; hold on; axis equal;
% xlabel('x'); ylabel('y');
% title('Figure 1: Configuration Sketch');
% 
% % Draw rectangle D
% rectangle('Position', [x1, y1, x2 - x1, y2 - y1], 'EdgeColor', 'b', 'LineWidth', 2);
% text(x2, -lambda/4, 'Domain D', 'HorizontalAlignment', 'center');
% 
% % Mark origin
% plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
% text(0, -lambda/5, '(0,0)', 'HorizontalAlignment', 'right');
% 
% % Mark right corner
% plot(lambda, lambda, 'ko', 'MarkerFaceColor', 'k');
% text(lambda+ lambda/2, lambda+lambda/4, '({\lambda},{\lambda})', 'HorizontalAlignment', 'right');
% 
% % Plot source location
% plot(sx, sy, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% text(sx, sy + lambda/4, 'Source s', 'HorizontalAlignment', 'center');
% 
% % Adjust axis
% %xlim([-lambda, sx + lambda]);
% ylim([-sy/10, sy + 2*lambda]);
% set(gca, 'YDir', 'reverse');  % because y increases downward
% grid on

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
u_inc0 = -1j/4 * besselh(0, 2, kb_0 * R);  % H_0^(2)(kb * R)

%% STEP 8
xc = lambda / 4;          
yc = lambda / 4;           
radius = lambda / 6;       % radius of the object
%side = lambda / 3;         % side length of the square object

% Define wavenumber field k(p)
k_object = 1.5* kb_0;       % Higher wavenumber inside the object
k_field = kb_0 * ones(Ny, Nx);  % Initialize as background

% Set k(p) > kb inside the circular object
object_mask = ((X - xc).^2 + (Y - yc).^2) <= radius^2;

%%Set k(p) > kb inside the square object
%object_mask = abs(X - xc) <= side/2 & abs(Y - yc) <= side/2;
%object_mask = abs(X - xc) <= side/2 & abs(Y - yc) <= lambda/40; %horizontal line
%object_mask = abs(X - xc) <= lambda/40 & abs(Y - yc) <= side/2; %vertical line

k_field(object_mask) = k_object;

% Define non-negative contrast function chi(ro)
chi = (k_field / kb_0).^2 - 1;

% STEP 9
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

% Plot source location
plot(sx, sy, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(sx, sy + lambda/4, 'Source s', 'HorizontalAlignment', 'center');

% Define the reciever domain 1
x_D_rec_start = -lambda; % x-coordinate of the start point of receiver domain
x_D_rec_end = 2*lambda; % x-coordinate of the end point of receiver domain
y_D_rec = 1.5*lambda; % y-coordinate of the receiver domain (does not change)
% Draw the D^{rec} in Figure 1
scatter(x_D_rec_start:h: x_D_rec_end, y_D_rec);
text(2.1*lambda, y_D_rec, 'Receiver Domain D^{rec}', 'HorizontalAlignment', 'left');

% Adjust axis
ylim([-1.5*lambda-sy/10, sy + 2*lambda]);
set(gca, 'YDir', 'reverse');  % because y increases downward
grid on
%% STEP 11
% Introduce grid on D^{rec}
%x_D_rec_grids = x_D_rec_start : 30*h : x_D_rec_end;
x_D_rec_grids = x_D_rec_start : 3*h : x_D_rec_end;

M = length(x_D_rec_grids);

%% STEP 13
N = size(chi,1) * size(chi,2);
chi_vec = reshape(chi, [N,1]); % contrast function chi(ro) vector
%% STEP 14
% Build the system matrix A (M x N)
% Define voxel area and midpoint values of the grids (voxels)
area = h^2; % area of a voxel via uniformly grid spacing
x_fine = 0:(h/2):lambda;      % 41 points: spacing = lambda/40
y_fine = 0:(h/2):lambda;
x_midpoint = x_fine(1:2:end);      % take every second value → 21 points
y_midpoint = y_fine(1:2:end);

[X_midpoint, Y_midpoint] = meshgrid(x_midpoint, y_midpoint);
kb_vals = linspace(1,10, 30);
A_all = zeros(M*length(kb_vals), N);

for idx = 1: length(kb_vals)
    kb = kb_vals(idx);
    Green_dist = zeros(M,N);
    for m = 1:M  
        % Compute |ro_m - ro_n'| at each grid point, ro_m -> receiver
        % locations, ro_n' -> midpoint locations of the voxels
        row_Green = sqrt((x_D_rec_grids(m) - X_midpoint).^2 + (y_D_rec - Y_midpoint).^2); 
        Green_dist(m,:) = reshape(row_Green, [1,N]);
    end
    
    % Compute the Green function (2D array) using Hankel function of the second kind, order 0
    Green_func = -1j/4 * besselh(0, 2, kb * Green_dist);  
    G_s = kb^2 * area * Green_func;
    % Compute the incident field (2D array) using Hankel function of the second kind, order 0
    u_inc = -1j/4 * besselh(0, 2, kb * R);  % H_0^(2)(kb * R)
    u_inc_vec = reshape(u_inc, [N,1]);
    A = G_s * diag(u_inc_vec);


    % % Compute the Green function (2D array) using Hankel function of the second kind, order 0
    % Green_func = -1j/4 * besselh(0, 2, kb * Green_dist);  
    % G_s = kb^2 * area * Green_func;
    % % Compute the incident field (2D array) using Hankel function of the second kind, order 0
    % u_inc = -1j/4 * besselh(0, 2, kb * R);  % H_0^(2)(kb * R)
    % % Create Kronecker product
    % K = kron(u_inc.', eye(Nx));  % u_inc.' is transpose (not conjugate transpose)
    % A = G_s * K;

    A_all(M*(idx-1)+1:M*idx, :) = A;
end

rank(A_all, 1e-3)
%% STEP 15

sing_vals = svd(A_all);

% Plot singular values
figure()
plot(1:length(sing_vals), sing_vals, 'x', 'MarkerSize', 10, 'LineWidth', 4); 
xlabel('Index', 'FontSize', 14);
ylabel('Singular Value', 'FontSize', 14);
title('Singular Values of A', 'FontSize', 16);
%set(gca, 'FontSize', 12);
set(gca,'YScale', 'log');
grid on;

% STEP 16
% Compute the scattered field
u_sc = A_all * chi_vec;

% STEP 17
% Find minimum norm solution to contrast function vector

%[U_r, S_r, V_r] = svd(A);
chi_vec_mn = pinv(A_all) * u_sc;

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

%% STEP 21
% Add noise
SNR_dB = 30; % e.g., 20 dB SNR
signal_power = norm(u_sc)^2 / length(u_sc);
noise_power = signal_power / (10^(SNR_dB/10));
noise = sqrt(noise_power) * randn(size(u_sc));  % Gaussian noise
u_sc_noisy = u_sc + noise;

chi_vec_mn_noisy = pinv(A_all) * u_sc_noisy;

chi_mn_noisy = reshape(real(chi_vec_mn_noisy), [Nx, Ny]);
chi_mn_noisy = max(chi_mn_noisy, 0);
% Plot
figure;
imagesc(x, y, chi_mn_noisy);
axis equal tight;
colorbar;
title('\chi_{mn}(\rho): Minimum Norm Contrast Function From Noisy Incident Field');
xlabel('x'); ylabel('y');
