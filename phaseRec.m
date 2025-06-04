function [holo_rec2, phaseImg, compSize] = phaseRec(filename, dx, dy, lambda, region, save)
            arguments
                filename {mustBeFile}
                dx {mustBeNumeric} 
                dy {mustBeNumeric} 
                lambda {mustBeNumeric}
                region {mustBeNumeric}
                save {mustBeNumeric}
            end
            
            disp(filename)
            disp(dx);
            disp(dy);
            disp(lambda);
            disp(region);
            disp(save);
            
            [path,name,ext] = fileparts(filename);

% Read hologram image and convert it to double precision
holo=double(imread(filename));
%imagesc(holo); % displays input image after coverted to doubles
% If the image is RGB, use only one channel (assume grayscale image)

% holo = holo(:,:,1);

% Get the size of the hologram
[N,M] = size(holo);
% Create a meshgrid for the hologram
[m,n] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);

% Calculate the Fourier Transform of the hologram and shift the zero-frequency component to the center
ft_holo = fftshift(fft2(fftshift(holo)));
figure, imagesc(log(abs(ft_holo).^2)), title 'Fourier tran' % idk why abs(log)^2
imagesc(log(abs(ft_holo).^2))
% Initialize a filter with zeros
filter = zeros(N,M);

% Create a filter mask for the desired region % make switch
if region==1
    filter(50:round(N/2-(N*0.19)),round(M/2+(M*0.18)):M-50) = 1; % 1nd quadrant
elseif region==2
    filter(1:round(N/2-(N*0.18)),1:round(M/2-(M*0.18))) = 1;  % 2nd quadrant
elseif region==3
    filter(round(N/2+(N*0.15)):N-40,40:round(M/2-(M*0.18))) = 1; % 3nd quadrant
else
    filter(round(N/2+(N*0.075)):N-500,round(M/2+(M*0.07)):M-600) = 1; % 4nd quadrant
end


% Apply the filter to the Fourier Transform of the hologram
ft_filtered_holo = ft_holo .* filter; %Max & min (from 10by10 conner)
filtered_spect = log(abs(ft_filtered_holo).^2);
figure, imagesc(filtered_spect), title 'Region filter'% find max?
mAx = 0; %max ftholo


% Threshold the phase image
BW = imbinarize(filtered_spect, 12);
figure;colormap gray;imagesc(BW);axis image, title 'BW'
[centers, radii, ~] = imfindcircles(BW,[10 800], ...
    ObjectPolarity="bright", Sensitivity=0.95);
[~,index] = max(radii);
centersStrong5 = centers(index,1:2);
radiiStrong5 = radii(index);
viscircles(centersStrong5, radiiStrong5,'EdgeColor','r');
[N,~] = size(filtered_spect);
resc = radiiStrong5; %radius to filter the object's spectrum. This is the radious of the +1/-1 terms
circ = ones(N);
for r = 1:N
    for p = 1:N
        if sqrt((r-centersStrong5(1,2))^2+(p-centersStrong5(1,1))^2)>resc
            circ(r,p) = 0;
        end
    end
end
filtrado = ft_filtered_holo.*circ;
figure;colormap gray;imagesc(log(abs(filtrado)));axis image, title 'filtrado circular'

filtered_spect = filtrado;
% Find the maximum value in the filtered spectrum
[~,idx] = max(filtered_spect(:));


% Define wavenumber
k = 2 * pi / lambda;

% Calculate the center frequencies for fx and fy
fx_0 = M/2;
fy_0 = N/2;

% Get the maximum values of fx and fy
[fy_max,fx_max] = ind2sub([N,M],idx);

% Define the step size for the search
step = 0.95;

% Initialize variables for the search
j = 0;

% Calculate the Inverse Fourier Transform of the filtered hologram
holo_rec = fftshift(ifft2(fftshift(filtered_spect)));

% Define the search range (G)
G = 3;
% Initialize flag for the search loop
fin = 0;

% Set initial values for fx and fy
fx = fx_max;
fy = fy_max;

% Initialize temporary search range
G_temp = G;

tic
% Loop to find the optimal fx and fy values
while fin == 0
  i = 0;
  j = j + 1;
  
  % Initialize the maximum sum (for thresholding)
  suma_maxima=0;
  
  % Nested loops for searching in the range of fx and fy
  for fy_tmp = fy-step*G_temp:step:fy+step*G_temp
    for  fx_tmp = fx-step*G_temp:step:fx+step*G_temp
      i = i + 1;
      
      % Calculate the metric for the current fx and fy
      suma = metric(holo_rec, fx_0, fy_0, fx_tmp, fy_tmp, lambda, M, N , dx, dy, m, n, k);
      
      % Update maximum sum and corresponding fx and fy if 
      % current sum is greater than the previous maximum
      if (suma > suma_maxima)
        x_max_out = fx_tmp;
        y_max_out = fy_tmp;
        suma_maxima=suma;
      end         
    end
  end
  
  % Update the temporary search range
  G_temp = G_temp - 1;
  
  % Check if the optimal values are found, set the flag to exit the loop
  if x_max_out == fx && y_max_out == fy
    fin = 1;        
  end
  
  % Update fx and fy for the next iteration
  fx = x_max_out;
  fy = y_max_out;  
end
toc
% Calculate the angles for the compensation wave
theta_x = asin((fx_0 - x_max_out) * lambda / (M * dx));
theta_y = asin((fy_0 - y_max_out) * lambda / (N * dy));

% Calculate the reference wave
ref = exp(1i * k * (sin(theta_x) * m * dx + sin(theta_y) * n * dy));

% Apply the reference wave to the hologram reconstruction
holo_rec2 = holo_rec .* (-ref);

compSize = size(holo_rec2);

% holo_rec2 = -holo_rec2;
Ampl = abs(holo_rec2);
phase = angle(holo_rec2);

phaseImg = phase;
ampImg = Ampl;

figure, imagesc(phase), colormap gray, title 'Phase image'
%figure, imagesc(Ampl), colormap gray

% Normalize the phase and convert it to uint8
phase = mat2gray(phase);
phase = uint8(phase * 255);

% Normalize the amplitude and convert it to uint8
Ampl = mat2gray(Ampl);
Ampl = uint8(Ampl * 255);

if save == true
    phaseName = strcat(path,name,"Phase_SHPC_output",ext);
    ampName = strcat(path,name,"Amplatude_SHPC_output",ext);

    imwrite(phase,gray(256),phaseName);
    imwrite(Ampl,gray(256),ampName);
end
end

function suma = metric(holo_rec, fx_0, fy_0, fx_tmp, fy_tmp, lambda, M, N , dx, dy, m, n, k)

   % Calculate the angles for the compensation wave
   theta_x = asin((fx_0 - fx_tmp) * lambda / (M * dx));
   theta_y = asin((fy_0 - fy_tmp) * lambda / (N * dy));

   % Calculate the reference wave
   ref = exp(1i * k * (sin(theta_x) * m * dx + sin(theta_y) * n * dy));

   % Apply the reference wave to the hologram reconstruction
   holo_rec2 = holo_rec .* -ref;
        
   % Calculate the phase of the hologram reconstruction
   phase = angle(holo_rec2);

   % Normalize the phase and convert it to uint8
   phase = mat2gray(phase);
   phase = uint8(phase * 255);

   % Threshold the phase image
   BW = imbinarize(phase, 0.13);

   % Calculate the sum of all elements in the resulting binary image
   suma = sum(sum(BW));
end


function [filtrado, circ] = findcircle(CF0_BW, CF0)
[centers, radii, ~] = imfindcircles(CF0_BW,[120 160], ObjectPolarity="bright", ...
    Sensitivity=0.95, EdgeThreshold=0.08);

[~, index] = max(radii);

centersStrong5 = centers(index,1:2);
radiiStrong5 = radii(index);

viscircles(centersStrong5, radiiStrong5,'EdgeColor','r');

[N,~] = size(CF0);

resc = radiiStrong5; %radius to filter the object's spectrum. This is the radious of the +1/-1 terms
circ = ones(N);
for r = 1:N
    for p = 1:N
        if sqrt((r-centersStrong5(1,2))^2+(p-centersStrong5(1,1))^2)>resc
            circ(r,p) = 0;
        end
    end
end

figure;colormap gray;imagesc(circ);axis image

filtrado = CF0.*circ;
end