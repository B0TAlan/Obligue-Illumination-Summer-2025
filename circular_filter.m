function [newFiltradoFT, circ0] = circular_filter(ComplexField, newCompleField_0, TH)
% TH = 15; %FOR BINARY
%TH = 12; %FOR CONTINOUS
CF0 = FT(ComplexField);
CF0_BW = imbinarize(log(abs(CF0).^2), TH);
% figure, imagesc(CF0_BW)
[filtrado0, circ0] = findcircle(CF0_BW, CF0);
newFiltradoFT = FT(newCompleField_0).*circ0;
% figure;colormap gray;imagesc(log(abs(newFiltradoFT.^2)));axis image
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