%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: functions_evaluation                                                  %
%                                                                              %
% The script contains all implemented function for SI_DHM_heuristic            %
%                                                                              %
% Authors: Raul Castaneda, Sofia Obando, Carlos Trujillo, Rene Restrepo,       %
%           Ana Doblas.                                                        %
% Applied Optics Group EAFIT univeristy                                        %
%                                                                              %
% Email: racastaneq@eafit.edu.co; adoblas@umassd.edu                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef function_heuristic
    methods(Static)

        function [holo,M,N,m,n] = holo_read(filename)

            % Function read the hologram
            % Inputs: filename  - file name with cmplete direction and
            %         extension
            % Output: holo - Matrix with the hologram
            %         M,N - dimsensions of the hologram
            %         m,n - space matrices of the size of holograms

            holo = double(imread(filename));
            holo = holo(:,:,1);
            [N,M] = size(holo);
            [n,m] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
        end


        function [holo_filter,holo_FT,fx_max,fy_max] = spatialFilter_SIDHM(holo,M,N,visual)

            % FROM OPHASE_R Function to apply a spatial filtering to the hologram
            % Inputs: holo - Matrix with the hologram
            %         M,N - dimsensions of the hologram
            %         m,n - space matrices of the size of holograms
            %         extension
            %         visual - logical value to visualize or not the
            %         filtering
            % Output: holo_filter - Filtered hologram
            %         holo_FT - Fourier Transform of the Hologram
            %         fx_max,fy_max - coordinates of the maximum within the
            %         +1 diffraction order
            fft_holo_1 = fftshift(fft2(fftshift(holo(:,:,1))));
            fft_holo_2 = fftshift(fft2(fftshift(holo(:,:,2))));

            % Initialize a filter with zeros
            filter = zeros(N,M);
            % Create a filter mask for the desired region
            filter(1:round(N/2-(N*0.18)),1:M) = 1; % 1nd quadrant

            ft_filtered_holo_1 = fft_holo_1 .* filter;
            filtered_spect1 = log(abs(ft_filtered_holo_1).^2);

            ft_filtered_holo_2 = fft_holo_2 .* filter;
            filtered_spect2 = log(abs(ft_filtered_holo_2).^2);

            BW = imbinarize(filtered_spect1, 12.5);
            BW = bwareafilt(BW, 1);
            BW_filled = imfill(BW, 'holes');
            edges = edge(BW_filled, 'Canny');

            [y,x] = find(edges);
            A = [-2*x, -2*y, ones(size(x))];
            b = -(x.^2 + y.^2);
            params = A\b;
            center(1) = params(1);
            center(2) = params(2);
            radius = ceil(sqrt(center(1)^2 + center(2)^2 - params(3)));
            % Resultados numéricos
            % fprintf('Centro detectado: (%.2f, %.2f)\n', center(1), center(2));
            % fprintf('Radio detectado: %.2f píxeles\n', radius);
            circ = ones(N);
            for r = 1:N
                for p = 1:N
                    if sqrt((r-center(2))^2+(p-center(1))^2)>radius
                        circ(r,p) = 0;
                    end
                end
            end
            Circ_filtered_spectrum1 = ft_filtered_holo_1.*circ;

            BW = imbinarize(filtered_spect2, 12);
            BW = bwareafilt(BW, 1);
            BW_filled = imfill(BW, 'holes');
            edges = edge(BW_filled, 'Canny');

            [y,x] = find(edges);
            A = [-2*x, -2*y, ones(size(x))];
            b = -(x.^2 + y.^2);
            params = A\b;
            center(1) = params(1);
            center(2) = params(2);
            radius = ceil(sqrt(center(1)^2 + center(2)^2 - params(3)));

            % Resultados numéricos
            % fprintf('Centro detectado: (%.2f, %.2f)\n', center(1), center(2));
            % fprintf('Radio detectado: %.2f píxeles\n', radius);
            circ2 = ones(N);
            for r = 1:N
                for p = 1:N
                    if sqrt((r-center(2))^2+(p-center(1))^2)>radius
                        circ2(r,p) = 0;
                    end
                end
            end
            Circ_filtered_spectrum2 = ft_filtered_holo_2.*circ2;

            holo_filter(:,:,1) = IFT(Circ_filtered_spectrum1);
            holo_filter(:,:,2) = IFT(Circ_filtered_spectrum2);

            holo_FT(:,:,1) = Circ_filtered_spectrum1;
            holo_FT(:,:,2) = Circ_filtered_spectrum2;

            % max values first peak
            maxValue_1 = max(max(abs(Circ_filtered_spectrum1)));
            [fy_max_1 fx_max_1] = find(abs(Circ_filtered_spectrum1) == maxValue_1);
            mask = ones(M,N);
            mask(fy_max_1 - 10:fy_max_1 + 10,fx_max_1 - 10:fx_max_1 + 10)= 0;
            fx_max_L = fx_max_1;
            fy_max_L = fy_max_1;
            Circ_filtered_spectrum1 = Circ_filtered_spectrum1 .* mask;

            maxValue_1 = max(max(abs(Circ_filtered_spectrum1)));
            [fy_max_1 fx_max_1] = find(abs(Circ_filtered_spectrum1) == maxValue_1);
            fx_max_D = fx_max_1(1);
            fy_max_D = fy_max_1(1);

            fy_max = [fx_max_L,fx_max_D];
            fx_max = [fy_max_L,fy_max_D];

            if visual == 'Yes'
                figure,imagesc(log(abs(fft_filter_holo).^2)),colormap(gray),title('FT Filter Hologram'),daspect([1 1 1])
                figure,imagesc((abs(holo_filter_1).^2)),colormap(gray),title('Filter Hologram'),daspect([1 1 1])
            end
        end


        function [cf] = costFunction_SIDHM(theta, FTHolo, fx_max,fy_max)

            % Function to apply The cost function to the first demodulation
            % Inputs: theta - angle for the demodulation
            %         FTHolo - Fourier Transform of the Hologram
            %         fx_max,fy_max - coordinates of the maximum within the
            %         +1 diffraction order
            % Output: cf - cost function for the first demodulation

            cf = 0;
            [M,N] = size(FTHolo);
            [Dtemp] = function_heuristic.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,1);
            cf =  abs(Dplus(fx_max(1),fy_max(1))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end


        function [cf] = costFunction_SIDHMII(theta, FTHolo, fx_max,fy_max)

            % Function to apply The cost function to the second demodulation
            % Inputs: theta - angle for the demodulation
            %         FTHolo - Fourier Transform of the Hologram
            %         fx_max,fy_max - coordinates of the maximum within the
            %         +1 diffraction order
            % Output: cf - cost function for the second demodulation

            cf = 0;
            [M,N] = size(FTHolo);
            [Dtemp] = function_heuristic.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,2);
            cf =  abs(Dplus(fx_max(2),fy_max(2))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end


        function [D] = demComp2SIDHM(theta, H)

            % Function to demodulate de periodic pattern of SI-DHM
            % Inputs: theta - pair of angles for demodulate
            %         H - Matrix with the two holograms, dimensions: (N,M,2)
            % Output: D - Metrix with the two demdulation G+ and G- dimensions: (N,M,2)

            [X, Y, no] = size(H);
            D = zeros(X,Y,no);
            M = 1/2*[exp(1i*0) exp(-1i*0);exp(1i*theta) exp(-1i*theta)];
            Minv = pinv(M);

            D(:,:,1) = Minv(1,1).*H(:,:,1) + Minv(1,2).*H(:,:,2);
            D(:,:,2) = Minv(2,1).*H(:,:,1) + Minv(2,2).*H(:,:,2);
        end


        function  [circ, ref] = phase_rec(filename, dx, dy, lambda)
            % Function to compensate the phase image
            % Inputs: Hologram - Hologram to compensate
            %         dx, dy  - pixel size in um
            %         lambda - wavelenght in um
            %         region - region of the cartessian plane where
            %         the +1  order is located
            %         save - parameter to save the results
            % Output: ref - numerical reference wave to compensate the hologram

            % Main function to perform phase reconstruction
            holo = filename;

            % Get the size of the hologram
            [N,M] = size(holo);
            % Create a meshgrid for the hologram
            [m,n] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);

            % Calculate the Fourier Transform of the hologram and shift the zero-frequency component to the center
            ft_holo = fftshift(fft2(fftshift(holo)));

            % Initialize a filter with zeros
            filter = zeros(N,M);
            % Create a filter mask for the desired region
            filter(2:round(N/2-(N*0.18)),2:M) = 1; % 1nd quadrant

            % Apply the filter to the Fourier Transform of the hologram
            ft_filtered_holo = ft_holo .* filter;
            filtered_spect = log(abs(ft_filtered_holo).^2);

            BW = imbinarize(filtered_spect, 10);
            BW = bwareafilt(BW, 1);
            BW_filled = imfill(BW, 'holes');
            edges = edge(BW_filled, 'Canny');

            [y,x] = find(edges);
            A = [-2*x, -2*y, ones(size(x))];
            b = -(x.^2 + y.^2);
            params = A\b;
            center(1) = params(1);
            center(2) = params(2);
            radius = ceil(sqrt(center(1)^2 + center(2)^2 - params(3)));
            % Resultados numéricos
            % fprintf('Centro detectado: (%.2f, %.2f)\n', center(1), center(2));
            % fprintf('Radio detectado: %.2f píxeles\n', radius);

            circ = ones(N);
            for r = 1:N
                for p = 1:N
                    if sqrt((r-center(2))^2+(p-center(1))^2)>radius
                        circ(r,p) = 0;
                    end
                end
            end
            Circ_filtered_spectrum1 = ft_filtered_holo.*circ;



            [~,idx] = max(Circ_filtered_spectrum1(:));
            % Get the maximum values of fx and fy
            [fy_max,fx_max] = ind2sub([N,M],idx);
            holo_rec = fftshift(ifft2(fftshift(Circ_filtered_spectrum1)));
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
            % figure, imagesc(phase), colormap gray, title 'Phase reconstruida'

            FT_holo = FT(holo_rec2);
            % figure, imagesc(log(abs(FT_holo).^2)), colormap gray, title 'FT centrado'


        end


        function suma = metric(holo_rec, fx_0, fy_0, fx_tmp, fy_tmp, lambda, M, N , dx, dy, m, n, k)

            % Function to evaluate the metric to find the best compensation possible
            % Inputs: holo_rec - filtered hologram
            %         fx_0, fy_0 - central frequency coordinates
            %         fx_tmp, fy_tmp - temporal frequency coordinates
            %         lambda - wavelenght in um
            %         M, N - dimsensions of the hologram
            %         dx, dy - pixel size in um
            %         m, n - space matrices of the size of holograms
            %         k - wave number
            % Output: suma - value of the metric

            % Calculate the angles for the compensation wave
            theta_x = asin((fx_0 - fx_tmp) * lambda / (M * dx));
            theta_y = asin((fy_0 - fy_tmp) * lambda / (N * dy));

            % Calculate the reference wave
            ref = exp(1i * k * (sin(theta_x) * m * dx + sin(theta_y) * n * dy));

            % Apply the reference wave to the hologram reconstruction
            holo_rec2 = holo_rec .* ref;

            % Calculate the phase of the hologram reconstruction
            phase = angle(holo_rec2);

            % Normalize the phase and convert it to uint8
            phase = mat2gray(phase);
            phase = uint8(phase * 255);

            % Threshold the phase image
            BW = imbinarize(phase, 0.1);

            % Calculate the sum of all elements in the resulting binary image
            suma = sum(sum(BW));
        end

    end
end
