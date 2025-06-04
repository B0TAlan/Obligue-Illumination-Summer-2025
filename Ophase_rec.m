function  holo_rec2 = phase_rec(filename, dx, dy, lambda, region)

arguments
                filename {mustBeFile}
                dx {mustBeNumeric} 
                dy {mustBeNumeric} 
                lambda {mustBeNumeric}
                region {mustBeNumeric}
            end

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
