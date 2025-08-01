function  [ref] = Sphase_rec(app,filename, dx, dy, lambda, circ, save)

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

            Circ_filtered_spectrum1 = ft_filtered_holo.*circ;

            [~,idx] = max(Circ_filtered_spectrum1(:));
            % Get the maximum values of fx and fy
            [fy_max,fx_max] = ind2sub([N,M],idx);

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
            if save == 1 % copy into flow (phase sum)
                fol = string(datetime("today","Format","uuuu-MM-dd"));
                met = strcat("dx:",string(dx),"\n","dy: ",string(dy),"\n","region: ",string(region),"\n","lambda: ",string(lambda));
                mkdir(fol);
                phaseName = strcat(fol,"\",name,"_OI_DHM_output",ext);
                tN = strcat(name,"_Params.txt");
                % make txt with params 
                writelines(met,tN,"WriteMode","overwrite");
                imwrite(phase,phaseName);
            end
        end

        
        function Ft_norm = normalization_FT(app,Given_FT)
            Ft_Gminus = fftshift(fft2(fftshift(Given_FT)));
            Ft_norm = (Ft_Gminus - min(min(Ft_Gminus)))/(max(max(Ft_Gminus)) - min(min(Ft_Gminus)));
        end
    end