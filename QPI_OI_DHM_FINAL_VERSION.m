% Clear workspace, close all figures, and clear command window
clear all
close all
clc

% Define file path and filename
path = 'C:\Users\sobandovasquez\OneDrive - University of Massachusetts Dartmouth\PhD\PUBLICACIONES\Cost-efective NanoLive\SIMULATIONS\';
ext = '.bmp';

% Define dx and dy (pixel pitch in um)
dx = 4.65; 
dy = 4.65;
% Define wavelength 
lambda = 0.532;
% Define the Cartessian region for the spatial filter 
region = 2;


%%
% 0 DEG StarBi_300_Thick_250
file = 'Star_0_Thick_250'
% Concatenate file path and extension
filename = strcat(path,file);
filename = strcat(filename,ext);
% Call phase_rec function to recover complex fiel of the sample
ComplexField0 = phase_rec(filename, dx, dy, lambda, region);
[N,M] = size(ComplexField0);
zero_phase = phase_normalization(ComplexField0);
newCompleField_0 = zero_phase;
% % 
phase_0deg = angle(newCompleField_0);
ns = 1.52; % borosilicato
nim = 1; % air
thickness_0 = (phase_0deg * lambda)/ (2* pi * (ns-nim));
% % 
[filtrado0, circ0] = circular_filter(ComplexField0, newCompleField_0);


%% Checking of the thickness

%profiles taken
xi = [M/2, M];
% yi = [N/2, N/2-60];
yi = [N/2, N/2];
perfil_maximo = improfile(thickness_0, xi, yi);
% xii = [512, 1021];
% yii = [512, 430];
% perfil_minimo = improfile(thickness_0, xii, yii);

%Check the profiles taken
figure(8), imagesc(thickness_0), colormap hot, title '0 deg thickness'
hold on
plot(xi, yi, 'c-', 'LineWidth', 2); % 'c-' for cyan line


%plot of the profiles
figure(6)
plot(perfil_maximo)



%% 60 DEG
file = 'Star_60_Thick_250'
% Concatenate file path and extension
filename = strcat(path,file);
filename = strcat(filename,ext);
% Call phase_rec function to recover complex fiel of the sample
ComplexField60 = phase_rec(filename, dx, dy, lambda, region);
zero_phase = phase_normalization(ComplexField60);
newCompleField_60 = zero_phase;
[filtrado60, circ60] = circular_filter(ComplexField60, newCompleField_60);
phase_60deg = angle(newCompleField_60);
thickness_60 = (phase_60deg * lambda)/ (2* pi * (ns-nim));
figure, imagesc(thickness_60), colormap hot, title '60 deg'

%% 120 DEG
file = 'Star_120_Thick_250'
% Concatenate file path and extension
filename = strcat(path,file);
filename = strcat(filename,ext);
% Call phase_rec function to recover complex fiel of the sample
ComplexField120 = phase_rec(filename, dx, dy, lambda, region);
zero_phase = phase_normalization(ComplexField120);
newCompleField_120 = zero_phase;
[filtrado120,circ120]=circular_filter(ComplexField120,newCompleField_120);
phase_120deg = angle(newCompleField_120);
thickness_120 = (phase_120deg * lambda)/ (2* pi * (ns-nim));
figure, imagesc(thickness_120), colormap hot, title '120 deg'

%% 180 DEG
file = 'Star_180_Thick_250'
% Concatenate file path and extension
filename = strcat(path,file);
filename = strcat(filename,ext);
% Call phase_rec function to recover complex fiel of the sample
ComplexField180 = phase_rec(filename, dx, dy, lambda, region);
zero_phase = phase_normalization(ComplexField180);
newCompleField_180 = zero_phase;
[filtrado180,circ180]=circular_filter(ComplexField180,newCompleField_180);
phase_180deg = angle(newCompleField_180);
thickness_180 = (phase_180deg * lambda)/ (2* pi * (ns-nim));
figure, imagesc(thickness_180), colormap hot, title '180 deg'

%% 240 DEG
file = 'Star_240_Thick_250'
% Concatenate file path and extension
filename = strcat(path,file);
filename = strcat(filename,ext);
% Call phase_rec function to recover complex fiel of the sample
ComplexField240 = phase_rec(filename, dx, dy, lambda, region);
zero_phase = phase_normalization(ComplexField240);
newCompleField_240 = zero_phase;
[filtrado240,circ240]=circular_filter(ComplexField240,newCompleField_240);
phase_240deg = angle(newCompleField_240);
thickness_240 = (phase_240deg * lambda)/ (2* pi * (ns-nim));
figure, imagesc(thickness_240), colormap hot, title '240 deg'

%% 300 DEG
file = 'Star_300_Thick_250'
% Concatenate file path and extension
filename = strcat(path,file);
filename = strcat(filename,ext);
% Call phase_rec function to recover complex fiel of the sample
ComplexField300 = phase_rec(filename, dx, dy, lambda, region);
zero_phase = phase_normalization(ComplexField300);
newCompleField_300 = zero_phase;
[filtrado300,circ300]=circular_filter(ComplexField300,newCompleField_300);
phase_300deg = angle(newCompleField_300);
thickness_300 = (phase_300deg * lambda)/ (2* pi * (ns-nim));
figure, imagesc(thickness_300), colormap hot, title '300 deg'

%% Composition of the flower
SUMA = filtrado0 + filtrado60 + filtrado120 + filtrado180 + filtrado240 + filtrado300;
figure, imagesc((abs(SUMA).^0.1)), colormap gray


%% Inverse of the flower
suma = IFT(SUMA);
phase_SRDHM = angle(suma);
figure, imagesc(phase_SRDHM), colormap hot, title 'flor sin normalizar en frecuencias'

%% Mask for the normalization
sum_circles = circ0 + circ60 + circ120 + circ180 + circ240 +circ300 ;
figure(20), imagesc(sum_circles), colormap gray

%% Inverse of the mask
% Calcular los pesos de normalización basados en las zonas de superposición
normalization_weights = 1 ./ sum_circles; % Inverso de la suma en cada punto
normalization_weights(isinf(normalization_weights)) = 0; % Manejar división por cero
figure(200), imagesc(normalization_weights), colormap gray
%
%% Normalizar la matriz sumada
normalized_circles = SUMA .* normalization_weights;
figure(202), imagesc(log(abs(normalized_circles).^2)), colormap gray


%% Normalized phase
suma = IFT(normalized_circles);
% suma = IFT(SUMA);
phase_SRDHM = angle(suma);
figure(21), imagesc(phase_SRDHM), colormap hot


thickness_total_SR = (phase_SRDHM * lambda)/ (2* pi * (ns-nim));
figure, imagesc(thickness_total_SR), colormap hot, title 'Thickness total SR'


%% Checking the thickness again

% profiles taken
xi = [M/2, M];
% yi = [N/2, N/2-70];

yi = [N/2+1, N/2];
perfil_maximo = improfile(thickness_total_SR, xi, yi);
% xii = [512, 1021];
% yii = [512, 430];
% perfil_minimo = improfile(thickness_0, xii, yii);

% Check the profiles taken
figure(8), imagesc(thickness_total_SR), colormap hot, title 'SR thickness'
hold on
plot(xi, yi, 'c-', 'LineWidth', 2); % 'c-' for cyan line
% hold on
% plot(xii, yii, 'g-', 'LineWidth', 2); % 'c-' for cyan line

%plot of the profiles
figure(6)
plot(perfil_maximo)
% hold on 
% plot((perfil_minimo))




%%
% DHM StarBi_300_Thick_250
file = 'Star_DHM_Thick_250'
% Concatenate file path and extension
filename = strcat(path,file);
filename = strcat(filename,ext);
% Call phase_rec function to recover complex fiel of the sample
ComplexFieldDHM = phase_rec(filename, dx, dy, lambda, region);
[N,M] = size(ComplexFieldDHM);
zero_phase = phase_normalization(ComplexFieldDHM);
newCompleField_DHM = zero_phase;
% % 
phase_DHMdeg = angle(newCompleField_DHM);
ns = 1.52; % borosilicato
nim = 1; % air
thickness_DHM = (phase_DHMdeg * lambda)/ (2* pi * (ns-nim));


% % 
[filtradoDHM, circDHM] = circular_filter(ComplexFieldDHM, ComplexFieldDHM);

%%
figure, imagesc(thickness_DHM)





%% Circular profiles
% 
radius_1 = 500; %px 
radius_2 = 200; %px
radius_3 = 20; %px 

% 
% radius_1 = 1500; %px 
% radius_2 = 500; %px
% radius_3 = 250; %px 

[theta_1, profile1] = perfil_circular(thickness_total_SR, radius_1);
[theta_2, profile2] = perfil_circular(thickness_total_SR, radius_2);
[theta_3, profile3] = perfil_circular(thickness_total_SR, radius_3);

[theta_DHM, profileDHM] = perfil_circular(thickness_DHM, radius_1);
[theta_DHM2, profileDHM2] = perfil_circular(thickness_DHM, radius_2);
[theta_DHM3, profileDHM3] = perfil_circular(thickness_DHM, radius_3);
%%
figure(20);
plot(rad2deg(theta_1), profile1, 'LineWidth', 1.5);
hold on
plot(rad2deg(theta_2), profile2, 'LineWidth', 1.5);
hold on 
plot(rad2deg(theta_3), profile3, 'LineWidth', 1.5);
legend('SR r = 500 px' , 'SR r = 200 px', 'SR r = 20 px')
title 'Perfiles circulares SR'


figure(21)
plot(rad2deg(theta_DHM), profileDHM, 'LineWidth', 1.5);
hold on
plot(rad2deg(theta_DHM2), profileDHM2, 'LineWidth', 1.5);
hold on
plot(rad2deg(theta_DHM3), profileDHM3, 'LineWidth', 1.5);
legend('DHM r = 500 px' , 'DHM r = 200 px', 'DHM r = 20 px')
title 'Perfiles circulares DHM'


%% functions


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



function [profile_norm, valid_idx, arc_length] = get_profile(img, radius)
% Get image size and define the center
[rows, cols] = size(img);
center = round([rows+1, cols+1] / 2)

% Define the inner and outer radius for the ring
inner_radius = radius;  % Adjust as needed
outer_radius = inner_radius +2; % Ring thickness = outer_radius - inner_radius
mean_radius = (inner_radius + outer_radius) / 2; % Approximate arc radius

% Define angular resolution (more points for smoother plot)
theta = linspace(0, 2*pi, 360); 

% Initialize intensity storage
profile = zeros(1, length(theta));

% Loop over multiple radii to get mean intensity inside the ring
for r = inner_radius:outer_radius
    % Compute coordinates for the current radius
    x_coords = round(center(2) + r * cos(theta));
    y_coords = round(center(1) + r * sin(theta));
    
    % Keep only valid coordinates inside the image
    valid_idx = (x_coords > 0 & x_coords <= cols) & (y_coords > 0 & y_coords <= rows);
    x_coords = x_coords(valid_idx);
    y_coords = y_coords(valid_idx);
    theta_valid = theta(valid_idx); % Keep angles that are valid
    
    % Extract and accumulate intensity values
    temp_profile = zeros(1, length(x_coords));
    for i = 1:length(x_coords)
        temp_profile(i) = img(y_coords(i), x_coords(i));
    end
    
    % Accumulate values for mean calculation
    profile(valid_idx) = profile(valid_idx) + temp_profile;
end

% % Display the original image with the ring overlaid
figure;
imshow(img, []);
hold on;
plot(center(2) + inner_radius * cos(theta), center(1) + inner_radius * sin(theta), 'b', 'LineWidth', 1);
plot(center(2) + outer_radius * cos(theta), center(1) + outer_radius * sin(theta), 'r', 'LineWidth', 1);
title('Image with Ring Overlay');

% Compute mean intensity in the ring
profile_norm = profile / (outer_radius - inner_radius +1);
% Convert angles to arc length in pixels
arc_length = mean_radius * theta_valid; % s = r * theta
end



function [newFiltradoFT, circ0] = circular_filter(ComplexField, newCompleField_0)
% TH = 15; %FOR BINARY
TH = 12; %FOR CONTINOUS
CF0 = FT(ComplexField);
CF0_BW = imbinarize(log(abs(CF0).^2), TH);
% figure, imagesc(CF0_BW)
[filtrado0, circ0] = findcircle(CF0_BW, CF0);
newFiltradoFT = FT(newCompleField_0).*circ0;
% figure;colormap gray;imagesc(log(abs(newFiltradoFT.^2)));axis image
end 

function ComplexFieldNorm = phase_normalization(ComplexField)
% phase = angle(ComplexField);
% min_phase_original = min(min(phase));
% relative_min = phase(1,1);
% zero_phase = phase - relative_min;
% mini_phase_2 = min(min(zero_phase));


phase = angle(ComplexField);
min_phase_original = min(min(phase));
ComplexFieldNorm = ComplexField.*exp(-1j*min_phase_original);
end

function [filtrado, circ] = findcircle(CF0_BW, CF0)
%Threshold the phase image
BW = CF0_BW;
[N,M] = size(CF0_BW);
figure;colormap gray;imagesc(BW);axis image, title 'BW'

% Limpiar la imagen, dejando solo la región más grande (círculo)
BW = bwareafilt(BW, 1);

% Rellenar el círculo binario para obtener un borde continuo y robusto
BW_filled = imfill(BW, 'holes');

% Detectar el borde del círculo rellenado
edges = edge(BW_filled, 'Canny');

% Extraer coordenadas del borde
[y,x] = find(edges);

% Ajuste robusto por mínimos cuadrados para círculo
A = [-2*x, -2*y, ones(size(x))];
b = -(x.^2 + y.^2);
params = A\b;

% Parámetros del círculo
center(1) = params(1);
center(2) = params(2);
radius = sqrt(center(1)^2 + center(2)^2 - params(3));

% Visualización clara de resultados
figure; imshow(BW); hold on;
viscircles([center(1), center(2)], radius, 'EdgeColor','b');
plot(center(1), center(2), 'r+', 'MarkerSize', 20, 'LineWidth', 2);

% Resultados numéricos
fprintf('Centro detectado: (%.2f, %.2f)\n', center(1), center(2));
fprintf('Radio detectado: %.2f píxeles\n', radius);

circ = ones(N);
for r = 1:N
    for p = 1:N
        if sqrt((r-center(2))^2+(p-center(1))^2)>radius
            circ(r,p) = 0;
        end
    end
end

filtrado = CF0.*circ;
end

% Main function to perform phase reconstruction
function  holo_rec2 = phase_rec(filename, dx, dy, lambda, region)

% Read hologram image and convert it to double precision
holo = double(imread(filename));
% Get the size of the hologram
[N,M] = size(holo);
% Create a meshgrid for the hologram
[m,n] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
% Calculate the Fourier Transform of the hologram 
% and shift the zero-frequency component to the center
ft_holo = fftshift(fft2(fftshift(holo)));
figure(15), imagesc(log(abs(ft_holo).^2)), title 'FT hologram'

% Initialize a filter with zeros
filter = zeros(N,M);
% Create a filter mask for the desired region
if region==1
    filter(50:round(N/2-(N*0.19)),round(M/2+(M*0.18)):M-50) = 1; % 1nd quadrant
elseif region==2
    filter(1:round(N/2-(N*0.185)),1:round(M/2-(M*0.185))) = 1;  % 2nd quadrant
elseif region==3
    filter(round(N/2+(N*0.15)):N-40,40:round(M/2-(M*0.18))) = 1; % 3nd quadrant
else
    filter(round(N/2+(N*0.075)):N-500,round(M/2+(M*0.07)):M-600) = 1; % 4nd quadrant
end
% Apply the filter to the Fourier Transform of the hologram
ft_filtered_holo = ft_holo .* filter;
filtered_spect = log(abs(ft_filtered_holo).^2);
figure(1), imagesc(filtered_spect)

%Threshold the phase image
% BW = imbinarize(filtered_spect, 15); %FOR BINARY
BW = imbinarize(filtered_spect, 12); %FOR CONTINOUS
figure;colormap gray;imagesc(BW);axis image, title 'BW'

% Limpiar la imagen, dejando solo la región más grande (círculo)
BW = bwareafilt(BW, 1);

% Rellenar el círculo binario para obtener un borde continuo y robusto
BW_filled = imfill(BW, 'holes');

% Detectar el borde del círculo rellenado
edges = edge(BW_filled, 'Canny');

% Extraer coordenadas del borde
[y,x] = find(edges);

% Ajuste robusto por mínimos cuadrados para círculo
A = [-2*x, -2*y, ones(size(x))];
b = -(x.^2 + y.^2);
params = A\b;

% Parámetros del círculo
center(1) = params(1);
center(2) = params(2);
radius = sqrt(center(1)^2 + center(2)^2 - params(3));

% Visualización clara de resultados
figure; imshow(BW); hold on;
viscircles([center(1), center(2)], radius, 'EdgeColor','b');
plot(center(1), center(2), 'r+', 'MarkerSize', 20, 'LineWidth', 2);

% Resultados numéricos
fprintf('Centro detectado: (%.2f, %.2f)\n', center(1), center(2));
fprintf('Radio detectado: %.2f píxeles\n', radius);

circ = ones(N);
for r = 1:N
    for p = 1:N
        if sqrt((r-center(2))^2+(p-center(1))^2)>radius
            circ(r,p) = 0;
        end
    end
end
filtered_spect = ft_filtered_holo.*circ;
figure;colormap gray;imagesc(log(abs(filtered_spect)));title 'filtrado circular'


[~,idx] = max(filtered_spect(:));
% Get the maximum values of fx and fy
[fy_max,fx_max] = ind2sub([N,M],idx);
holo_rec = fftshift(ifft2(fftshift(filtered_spect)));
% Find the maximum value in the filtered spectrum
% Define wavenumber
k = 2 * pi / lambda;
% Calculate the center frequencies for fx and fy
fx_0 = M/2;
fy_0 = N/2;

% Calculate the angles for the compensation wave
theta_x = asin((fx_0 - fx_max +1) * lambda / (M * dx));
theta_y = asin((fy_0 - fy_max +1) * lambda / (N * dy));
% Calculate the reference wave
ref = exp(1i * k * (sin(theta_x) * m * dx + sin(theta_y) * n * dy));
% Apply the reference wave to the hologram reconstruction
holo_rec2 = holo_rec .* (-ref);

phase = angle(holo_rec2);
figure, imagesc(phase), colormap gray, title 'Phase reconstruida'

FT_holo = FT(holo_rec2);
figure, imagesc(log(abs(FT_holo).^2)), colormap gray, title 'FT centrado'

end

