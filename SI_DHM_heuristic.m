%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: SI_DHM_heuristics                                                     %
%                                                                              %
% The script cperforms the demodulation and composition of SR images from      %
% SI-DHM hologram witha random phase shift between them                        %
%                                                                              %                                       
% Authors: Raul Castaneda, Sofia Obando, Carlos Trujillo, Rene Restrepo,       %
%           Ana Doblas.                                                        %
% Applied Optics Group EAFIT univeristy, Colombia                              % 
% and Optica Imaging Research Laboratory University of MAssachusetts Dartmouth % 
%                                                                              %
% Email: racastaneq@eafit.edu.co; adoblas@umassd.edu                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc% clean
close all% close all windows
clear all% clear of memory all variable

%% Read recorded holograms

% Lines to add folders for reading images and/or functions
% EAFIT computer
addpath('C:\Users\ASUS\OneDrive - University of Massachusetts Dartmouth\PhD\CONFERENCES\SPIE OPTICS AND PHOTONICS\ALAN POSTER\MATLAB\')
% addpath('C:\Users\racastaneq\Documents\MEGA\MEGAsync\RACQ\Universities\05 EAFIT\Research projects\2023\SI-SOFI\GitHub versions\Matlab\')

% Load the hologram
name1 =  'struct IMG\holo1_angle0.bmp'; 
name2 =  'struct IMG\holo2_angle0.bmp'; 

[h1,M,N,m,n] = function_heuristic.holo_read(name1);
[h2,M,N,m,n] = function_heuristic.holo_read(name2);

% Input parameters to record the holograms
lambda = 0.532;   % Source's wavelength in microns
k = 2*pi/lambda;  % Wavenumber in microns
dxy = 4.65;       % Pitch in X and Y direction in microns
fx_0 = M/2;
fy_0 = N/2;

% Fourier transform of the holograms
H1 = fftshift(fft2(fftshift(h1)));
H2 = fftshift(fft2(fftshift(h2)));

% Display the Fourier spectrum of the first hologram
figure; colormap gray; imagesc(abs(H1).^0.1); axis image;
%%
% Lines to implement the spatial filter using a circular mask
holo(:,:,1) = h1;
holo(:,:,2) = h2;
[holo_filter,holo_FT,fx_max,fy_max] = function_heuristic.spatialFilter_SIDHM(holo,M,N,'yes');


% minimization particleswarm
repeat_minimization = true;
while repeat_minimization
    % Minimization particleswarm
    % Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
    test_theta = randi([0 360])*pi/180;
    lb = 0;
    ub = 360*pi/180;

    options = optimoptions('particleswarm', 'Display', 'off', 'SwarmSize', 4, ...
        'MaxIterations', 500, 'FunctionTolerance', 1e-8);
    cf_particleswarm = @(test_theta) function_heuristic.costFunction_SIDHM(test_theta, holo_FT, fx_max, fy_max);

    tic
    [theta, cf_particleswarm_val] = particleswarm(@(params) function_heuristic.costFunction_SIDHM(params, holo_FT, fx_max, fy_max), 1, lb, ub, options);
    toc

    fprintf('theta 1: %f\n', theta);

    [Gdemod] = function_heuristic.demComp2SIDHM(theta, holo_FT);
    Gplus_demod = Gdemod(:,: ,1);
    figure, imagesc(log(abs(Gplus_demod).^2))
 
    % Ask the user if they want to repeat
    answer = input('Do you want to repeat the minimization? (y/n): ', 's');
    if lower(answer) ~= 'y'
        repeat_minimization = false;
    end
end
close all

gplus = fftshift(ifft2(fftshift(Gplus_demod)));
gplus = gplus * -1;
[circ, ref_gplus] = function_heuristic.phase_rec(gplus, dxy, dxy, lambda);
gplus1 = gplus .* ref_gplus;
Gplus1 = FT(gplus1);

% Find max in A_before
[~, idx_before] = max(log(abs(Gplus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gplus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gplus1(:)).^2)));
[row_after, col_after] = ind2sub(size(Gplus1), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ1 = circshift(circ, [delta_row, delta_col]);

FT_centrada = FT(gplus1).*circ1;
Centered_complex_field = IFT(FT_centrada);
zero_phase = (Centered_complex_field);

figure, imagesc((abs(FT(zero_phase)).^.2)), colormap gray

% Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
repeat_minimization = true;

while repeat_minimization
    %%% Minimization particleswarm
    % Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
    test_theta = randi([0 360])*pi/180;
    cf_particleswarm = function_heuristic.costFunction_SIDHMII(test_theta,holo_FT,fx_max,fy_max);

    tic
    [theta, cf_particleswarm] = particleswarm(@(params) function_heuristic.costFunction_SIDHMII(params,holo_FT,fx_max,fy_max), 1, lb, ub, options);
    toc

    fprintf('theta 2: %f\n', theta);

    [Gdemod] =  function_heuristic.demComp2SIDHM(theta, holo_FT);
    Gminus_demod = Gdemod(:,:,2);
    figure, imagesc(log(abs(Gminus_demod).^2))

    % Ask the user if they want to repeat
    answer = input('Do you want to repeat the minimization? (y/n): ', 's');
    if lower(answer) ~= 'y'
        repeat_minimization = false;
    end
end
close all

gminus = fftshift(ifft2(fftshift(Gminus_demod)));
[circ, ref_minus] = function_heuristic.phase_rec(gminus, dxy, dxy, lambda);
gminus1 = gminus .* ref_minus;
Gminus1 = FT(gminus1);
% Find max in A_before
[~, idx_before] = max(log(abs(Gminus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gminus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gminus1(:)).^2)));
[row_after, col_after] = ind2sub(size(Gminus1), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ2 = circshift(circ, [delta_row, delta_col]);

FT_centrada2 = FT(gminus1).*circ2;
Centered_complex_field2 = IFT(FT_centrada2);
zero_phase2 = (Centered_complex_field2);

% figure, imagesc(angle(zero_phase2))
figure, imagesc((abs(FT(zero_phase2)).^.2)), colormap gray
%%

% Load the hologram
name1 =  'struct IMG\holo1_angle60.bmp'; 
name2 =  'struct IMG\holo2_angle60.bmp'; 

[h1,M,N,m,n] = function_heuristic.holo_read(name1);
[h2,M,N,m,n] = function_heuristic.holo_read(name2);

% Fourier transform of the holograms
H1 = fftshift(fft2(fftshift(h1)));
H2 = fftshift(fft2(fftshift(h2)));

% Display the Fourier spectrum of the first hologram
figure; colormap gray; imagesc(abs(H1).^0.1); axis image;

% Lines to implement the spatial filter using a circular mask
holo(:,:,1) = h1;
holo(:,:,2) = h2;
[holo_filter,holo_FT,fx_max,fy_max] = function_heuristic.spatialFilter_SIDHM(holo,M,N,'Not');


% minimization particleswarm
repeat_minimization = true;
while repeat_minimization
    % Minimization particleswarm
    % Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
    test_theta = randi([0 360])*pi/180;
    lb = 0;
    ub = 360*pi/180;

    options = optimoptions('particleswarm', 'Display', 'off', 'SwarmSize', 4, ...
        'MaxIterations', 500, 'FunctionTolerance', 1e-8);
    cf_particleswarm = @(test_theta) function_heuristic.costFunction_SIDHM(test_theta, holo_FT, fx_max, fy_max);

    tic
    [theta, cf_particleswarm_val] = particleswarm(@(params) function_heuristic.costFunction_SIDHM(params, holo_FT, fx_max, fy_max), 1, lb, ub, options);
    toc

    fprintf('theta 1: %f\n', theta);

    [Gdemod] = function_heuristic.demComp2SIDHM(theta, holo_FT);
    Gplus_demod = Gdemod(:,: ,1);
    figure, imagesc(log(abs(Gplus_demod).^2))
 
    % Ask the user if they want to repeat
    answer = input('Do you want to repeat the minimization? (y/n): ', 's');
    if lower(answer) ~= 'y'
        repeat_minimization = false;
    end
end
close all

gplus = fftshift(ifft2(fftshift(Gplus_demod)));
gplus = gplus * -1;
[circ, ref_gplus] = function_heuristic.phase_rec(gplus, dxy, dxy, lambda);
gplus3 = gplus .* ref_gplus;
Gplus3 = FT(gplus3);

% Find max in A_before
[~, idx_before] = max(log(abs(Gplus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gplus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gplus3(:)).^2)));
[row_after, col_after] = ind2sub(size(Gplus3), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ3 = circshift(circ, [delta_row, delta_col]);

FT_centrada3 = FT(gplus3).*circ3;
Centered_complex_field3 = IFT(FT_centrada3);
zero_phase3 = (Centered_complex_field3);

figure, imagesc((abs(FT(zero_phase3)).^.2)), colormap gray

% Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
repeat_minimization = true;

while repeat_minimization
    %%% Minimization particleswarm
    % Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
    test_theta = randi([0 360])*pi/180;
    cf_particleswarm = function_heuristic.costFunction_SIDHMII(test_theta,holo_FT,fx_max,fy_max);

    tic
    [theta, cf_particleswarm] = particleswarm(@(params) function_heuristic.costFunction_SIDHMII(params,holo_FT,fx_max,fy_max), 1, lb, ub, options);
    toc

    fprintf('theta 2: %f\n', theta);

    [Gdemod] =  function_heuristic.demComp2SIDHM(theta, holo_FT);
    Gminus_demod = Gdemod(:,:,2);
    figure, imagesc(log(abs(Gminus_demod).^2))

    % Ask the user if they want to repeat
    answer = input('Do you want to repeat the minimization? (y/n): ', 's');
    if lower(answer) ~= 'y'
        repeat_minimization = false;
    end
end
close all

gminus = fftshift(ifft2(fftshift(Gminus_demod)));
[circ, ref_minus] = function_heuristic.phase_rec(gminus, dxy, dxy, lambda);
gminus1 = gminus .* ref_minus;
Gminus1 = FT(gminus1);
% Find max in A_before
[~, idx_before] = max(log(abs(Gminus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gminus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gminus1(:)).^2)));
[row_after, col_after] = ind2sub(size(Gminus1), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ4 = circshift(circ, [delta_row, delta_col]);

FT_centrada4 = FT(gminus1).*circ4;
Centered_complex_field4 = IFT(FT_centrada4);
zero_phase4 = (Centered_complex_field4);

% figure, imagesc(angle(zero_phase2))
figure, imagesc((abs(FT(zero_phase4)).^.2)), colormap gray


%%
% Load the hologram
name1 =  'struct IMG\holo1_angle120.bmp'; 
name2 =  'struct IMG\holo2_angle120.bmp'; 

[h1,M,N,m,n] = function_heuristic.holo_read(name1);
[h2,M,N,m,n] = function_heuristic.holo_read(name2);

% Fourier transform of the holograms
H1 = fftshift(fft2(fftshift(h1)));
H2 = fftshift(fft2(fftshift(h2)));

% Display the Fourier spectrum of the first hologram
figure; colormap gray; imagesc(abs(H1).^0.1); axis image;

% Lines to implement the spatial filter using a circular mask
holo(:,:,1) = h1;
holo(:,:,2) = h2;
[holo_filter,holo_FT,fx_max,fy_max] = function_heuristic.spatialFilter_SIDHM(holo,M,N,'Not');



% minimization particleswarm
repeat_minimization = true;
while repeat_minimization
    % Minimization particleswarm
    % Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
    test_theta = randi([0 360])*pi/180;
    lb = 0;
    ub = 360*pi/180;

    options = optimoptions('particleswarm', 'Display', 'off', 'SwarmSize', 4, ...
        'MaxIterations', 500, 'FunctionTolerance', 1e-8);
    cf_particleswarm = @(test_theta) function_heuristic.costFunction_SIDHM(test_theta, holo_FT, fx_max, fy_max);

    tic
    [theta, cf_particleswarm_val] = particleswarm(@(params) function_heuristic.costFunction_SIDHM(params, holo_FT, fx_max, fy_max), 1, lb, ub, options);
    toc

    fprintf('theta 1: %f\n', theta);

    [Gdemod] = function_heuristic.demComp2SIDHM(theta, holo_FT);
    Gplus_demod = Gdemod(:,: ,1);
    figure, imagesc(log(abs(Gplus_demod).^2))
 
    % Ask the user if they want to repeat
    answer = input('Do you want to repeat the minimization? (y/n): ', 's');
    if lower(answer) ~= 'y'
        repeat_minimization = false;
    end
end
close all

gplus = fftshift(ifft2(fftshift(Gplus_demod)));
gplus = gplus * -1;
[circ, ref_gplus] = function_heuristic.phase_rec(gplus, dxy, dxy, lambda);
gplus5 = gplus .* ref_gplus;
Gplus5 = FT(gplus5);

% Find max in A_before
[~, idx_before] = max(log(abs(Gplus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gplus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gplus5(:)).^2)));
[row_after, col_after] = ind2sub(size(Gplus5), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ5 = circshift(circ, [delta_row, delta_col]);

FT_centrada5 = FT(gplus5).*circ5;
Centered_complex_field5 = IFT(FT_centrada5);
zero_phase5 = (Centered_complex_field5);

figure, imagesc((abs(FT(zero_phase5)).^.2)), colormap gray

% Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
repeat_minimization = true;

while repeat_minimization
    %%% Minimization particleswarm
    % Blind demodulation to recover the two shifted object spectrum individually without prior knowledge of the phase in the SI pattern
    test_theta = randi([0 360])*pi/180;
    cf_particleswarm = function_heuristic.costFunction_SIDHMII(test_theta,holo_FT,fx_max,fy_max);

    tic
    [theta, cf_particleswarm] = particleswarm(@(params) function_heuristic.costFunction_SIDHMII(params,holo_FT,fx_max,fy_max), 1, lb, ub, options);
    toc

    fprintf('theta 2: %f\n', theta);

    [Gdemod] =  function_heuristic.demComp2SIDHM(theta, holo_FT);
    Gminus_demod = Gdemod(:,:,2);
    figure, imagesc(log(abs(Gminus_demod).^2))

    % Ask the user if they want to repeat
    answer = input('Do you want to repeat the minimization? (y/n): ', 's');
    if lower(answer) ~= 'y'
        repeat_minimization = false;
    end
end
close all

gminus = fftshift(ifft2(fftshift(Gminus_demod)));
[circ, ref_minus] = function_heuristic.phase_rec(gminus, dxy, dxy, lambda);
gminus6 = gminus .* ref_minus;
Gminus1 = FT(gminus6);
% Find max in A_before
[~, idx_before] = max(log(abs(Gminus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gminus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gminus1(:)).^2)));
[row_after, col_after] = ind2sub(size(Gminus1), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ6 = circshift(circ, [delta_row, delta_col]);

FT_centrada6 = FT(gminus6).*circ6;
Centered_complex_field6 = IFT(FT_centrada6);
zero_phase6 = (Centered_complex_field6);

% figure, imagesc(angle(zero_phase2))
figure, imagesc((abs(FT(zero_phase6)).^.2)), colormap gray



%%
Ft_norm1 = normalization_FT(FT(zero_phase));
Ft_norm2 = normalization_FT(FT(zero_phase2));
Ft_norm3 = normalization_FT(FT(zero_phase3));
Ft_norm4 = normalization_FT(FT(zero_phase4));
Ft_norm5 = normalization_FT(FT(zero_phase5));
Ft_norm6 = normalization_FT(FT(zero_phase6));

Gsim = Ft_norm1 + Ft_norm2+ Ft_norm3 + Ft_norm4 + Ft_norm5 + Ft_norm6;
% Gsim = FT(zero_phase) + FT(zero_phase2) + FT(zero_phase3) + FT(zero_phase4) + FT(zero_phase5) + FT(zero_phase6);
gsim = fftshift(ifft2(fftshift(Gsim)));

figure,imagesc(angle(gsim)),colormap(gray),colorbar,title('Phase sim'),daspect([1 1 1])

figure,imagesc((abs(Gsim).^.1)),colormap(gray),colorbar,title('FT SIM'),daspect([1 1 1])

%% Mask for the normalization
sum_circles = circ1 + circ2 +circ3 + circ4  + circ5 + circ6 ;
figure(20), imagesc(sum_circles), colormap gray

% Inverse of the mask
% Calcular los pesos de normalización basados en las zonas de superposición
normalization_weights = 1 ./ sum_circles; % Inverso de la suma en cada punto
normalization_weights(isinf(normalization_weights)) = 0; % Manejar división por cero
figure(200), imagesc(normalization_weights), colormap gray



%% Normalizar la matriz sumada
normalized_circles = Gsim .* normalization_weights;
figure(202), imagesc((abs(normalized_circles).^.1)), colormap gray

%% Normalized phase
suma = IFT(normalized_circles);
% suma = IFT(SUMA);
phase_SRDHM = angle(suma);
figure(21), imagesc(phase_SRDHM), colormap hot

ns = 1.52; % borosilicato
nim = 1; % air
thickness_total_SR = (phase_SRDHM * lambda)/ (2* pi * (ns-nim));

thickness_total_SR = thickness_total_SR - min(min(thickness_total_SR));
figure, imagesc(thickness_total_SR), colormap hot, title 'Thickness total SR'


imgNorm = mat2gray(thickness_total_SR);
[I,~] = gray2ind(imgNorm,256);   
cmap   = hot(256);              
imwrite(I, cmap, 'SI_DHM_phase_map_normal.bmp'); 

%% Circular profiles
radius_1 = 500; %px 
radius_2 = 200; %px
radius_3 = 30; %px 

[theta_1, profile1] = perfil_circular(thickness_total_SR, radius_1);
[theta_2, profile2] = perfil_circular(thickness_total_SR, radius_2);
[theta_3, profile3] = perfil_circular(thickness_total_SR, radius_3);

%%
figure(20);
plot(rad2deg(theta_1), profile1, 'LineWidth', 1.5);
hold on
plot(rad2deg(theta_2), profile2, 'LineWidth', 1.5);
hold on 
plot(rad2deg(theta_3), profile3, 'LineWidth', 1.5);
legend('SR r = 500 px' , 'SR r = 200 px', 'SR r = 30 px')
title 'Perfiles circulares SR'


%%



% ------------------------------------------------------------------------
%  FIGURA 1 – Perfiles SR
% ------------------------------------------------------------------------
figure('Name','Perfiles circulares – SR','Color','w');
tSR = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile
plot(theta_1, profile1,'LineWidth',1.2);
title(sprintf('SR – radius = %d px',radius_1));
xlim([0 2*pi])
ylabel('Height (um)');
grid on

nexttile
plot(theta_2, profile2,'LineWidth',1.2);
title(sprintf('SR – radius = %d px',radius_2));
xlim([0 2*pi])
ylabel('Height (um)');
grid on

nexttile
plot(theta_3, profile3,'LineWidth',1.2);
title(sprintf('SR – radius = %d px',radius_3));
xlim([0 2*pi])
ylabel('Height (um)');
xlabel('\theta (rad)');
grid on

%% functions 

function ComplexFieldNorm = phase_normalization(ComplexField)
phase = angle(ComplexField);
min_phase_original = min(min(phase));
ComplexFieldNorm = ComplexField.*exp(-1j*min_phase_original);
end


function [theta , profile] = perfil_circular(phase, circle_radius)

% Obtener dimensiones de la imagen
[rows, cols] = size(phase);
center = ([rows+1, cols+1] / 2); % Centro de la imagen

% Definir el radio del círculo a analizar (ajustar según necesidad)


% Definir los ángulos para el círculo
theta = linspace(0, 2*pi, 360); % 360 puntos para cubrir el círculo

% Obtener coordenadas del círculo
x_coords = round(center(2) + circle_radius * cos(theta));
y_coords = round(center(1) + circle_radius * sin(theta));

% Evitar índices fuera de la imagen
valid_idx = (x_coords > 0 & x_coords <= cols) & (y_coords > 0 & y_coords <= rows);
x_coords = x_coords(valid_idx);
y_coords = y_coords(valid_idx);
theta = theta(valid_idx); % Mantener los ángulos válidos

% Obtener valores de intensidad a lo largo del círculo
profile = zeros(size(theta));
for i = 1:length(theta)
    profile(i) = phase(y_coords(i), x_coords(i));
end

% % Display the original image with the ring overlaid
figure;
imshow(phase, []);
hold on;
plot(center(2) + circle_radius * cos(theta), center(1) + circle_radius * sin(theta), 'b', 'LineWidth', 1);
title('Image with Ring Overlay');
end


function Ft_norm = normalization_FT(Given_FT)
Ft_Gminus = Given_FT;
Ft_norm = (Ft_Gminus - min(min(Ft_Gminus)))/(max(max(Ft_Gminus)) - min(min(Ft_Gminus)));
end 