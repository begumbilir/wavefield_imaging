%% STEP 7:
% Draw object domain and the locate the source 
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

%%
h = lambda / 20; % Grid spacing
x = 0:h:lambda; % x coordinates (21 points)
y = 0:h:lambda; % y coordinates (21 points)
Nx = length(x);
Ny = length(y);

[X, Y] = meshgrid(x, y);   % Y increases downwards by default in images

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

%% Change the locaton of the source (move closer to obj domain D)
% Location of the source
sx = 1/2 * lambda;
sy =  2*lambda;

% Figure
figure; hold on; axis equal;
xlabel('x'); ylabel('y');
title('Figure 1: Configuration Sketch');

% Draw rectangle D
rectangle('Position', [x1, y1, x2 - x1, y2 - y1], 'EdgeColor', 'b', 'LineWidth', 2);
text(x2/2, -lambda/5, 'Domain D', 'HorizontalAlignment', 'center');

% Mark origin
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
text(0, -lambda/10, '(0,0)', 'HorizontalAlignment', 'right');

% Mark right corner
plot(lambda, lambda, 'ko', 'MarkerFaceColor', 'k');
text(lambda+ lambda/5, lambda+lambda/10, '({\lambda},{\lambda})', 'HorizontalAlignment', 'right');

% Plot source location
plot(sx, sy, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(sx, sy + lambda/4, 'Source s', 'HorizontalAlignment', 'center');
ylim([-sy/5, sy + 2*lambda]);
set(gca, 'YDir', 'reverse');  % because y increases downward
grid on

% Compute |ro - ro_s| at each grid point
R = sqrt((X - sx).^2 + (Y - sy).^2);

% Compute the incident field (2D array) using Hankel function of the second kind, order 0
u_inc_new = -1j/4 * besselh(0, 2, kb * R);  % H_0^(2)(kb * R)

% Plot real, imag and absolute values

figure;

% Real part
subplot(1,3,1);
imagesc(x, y, real(u_inc_new));
axis equal tight;
title('Re(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Imaginary part
subplot(1,3,2);
imagesc(x, y, imag(u_inc_new));
axis equal tight;
title('Im(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Magnitude
subplot(1,3,3);
imagesc(x, y, abs(u_inc_new));
axis equal tight;
title('|u^{inc}|');
xlabel('x'); ylabel('y');
colorbar;

% Reverse y-axis to match coordinate system orientation (y increases downward)
set(gca,'YDir','normal');


%% Change kb value to observe the outcomes
kb = 2; 
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

%%
h = lambda / 20; % Grid spacing
x = 0:h:lambda; % x coordinates (21 points)
y = 0:h:lambda; % y coordinates (21 points)
Nx = length(x);
Ny = length(y);

[X, Y] = meshgrid(x, y);   % Y increases downwards by default in images

% Compute |ro - ro_s| at each grid point
R = sqrt((X - sx).^2 + (Y - sy).^2);

% Compute the incident field (2D array) using Hankel function of the second kind, order 0
u_inc_2 = -1j/4 * besselh(0, 2, kb * R);  % H_0^(2)(kb * R)

% Plot real, imag and absolute values

figure;

% Real part
subplot(1,3,1);
imagesc(x, y, real(u_inc_2));
axis equal tight;
title('Re(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Imaginary part
subplot(1,3,2);
imagesc(x, y, imag(u_inc_2));
axis equal tight;
title('Im(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Magnitude
subplot(1,3,3);
imagesc(x, y, abs(u_inc_2));
axis equal tight;
title('|u^{inc}|');
xlabel('x'); ylabel('y');
colorbar;

% Reverse y-axis to match coordinate system orientation (y increases downward)
set(gca,'YDir','normal');

%% Change the locaton of the source (move closer to obj domain D)
% Location of the source
sx = 1/2 * lambda;
sy = 2*lambda;

% Figure
figure; hold on; axis equal;
xlabel('x'); ylabel('y');
title('Figure 1: Configuration Sketch');

% Draw rectangle D
rectangle('Position', [x1, y1, x2 - x1, y2 - y1], 'EdgeColor', 'b', 'LineWidth', 2);
text(x2/2, -lambda/5, 'Domain D', 'HorizontalAlignment', 'center');

% Mark origin
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
text(0, -lambda/10, '(0,0)', 'HorizontalAlignment', 'right');

% Mark right corner
plot(lambda, lambda, 'ko', 'MarkerFaceColor', 'k');
text(lambda+ lambda/5, lambda+lambda/10, '({\lambda},{\lambda})', 'HorizontalAlignment', 'right');

% Plot source location
plot(sx, sy, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(sx, sy + lambda/4, 'Source s', 'HorizontalAlignment', 'center');
ylim([-sy/5, sy + 2*lambda]);
set(gca, 'YDir', 'reverse');  % because y increases downward
grid on

% Compute |ro - ro_s| at each grid point
R = sqrt((X - sx).^2 + (Y - sy).^2);

% Compute the incident field (2D array) using Hankel function of the second kind, order 0
u_inc_new_2 = -1j/4 * besselh(0, 2, kb * R);  % H_0^(2)(kb * R)

% Plot real, imag and absolute values

figure;

% Real part
subplot(1,3,1);
imagesc(x, y, real(u_inc_new_2));
axis equal tight;
title('Re(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Imaginary part
subplot(1,3,2);
imagesc(x, y, imag(u_inc_new_2));
axis equal tight;
title('Im(u^{inc})');
xlabel('x'); ylabel('y');
colorbar;

% Magnitude
subplot(1,3,3);
imagesc(x, y, abs(u_inc_new_2));
axis equal tight;
title('|u^{inc}|');
xlabel('x'); ylabel('y');
colorbar;

% Reverse y-axis to match coordinate system orientation (y increases downward)
set(gca,'YDir','normal');