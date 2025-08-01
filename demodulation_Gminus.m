function [circ1,zero_phase6] = demodulation_Gminus(std_temp_gminus, angles,FTHolo, dxy, lambda, circ2)
posicion = find(std_temp_gminus == min(std_temp_gminus));
theta_final = [0, angles(posicion)];
[Gdemod] =  demComp2SIDHM(theta_final, FTHolo);
Gminus_demod = Gdemod(:,:,2);

gminus = fftshift(ifft2(fftshift(Gminus_demod)));
gminus = gminus*-1;
% ref_minus = phase_rec(gminus, dxy, dxy, lambda, region, false);
[ref_minus] = phase_rec(gminus, dxy, dxy, lambda, circ2);
Gminus = gminus .*ref_minus;

Gminus1 = FT(Gminus);
% Find max in A_before
[~, idx_before] = max(log(abs(Gminus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gminus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gminus1(:)).^2)));
[row_after, col_after] = ind2sub(size(Gminus1), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ1 = circshift(circ2, [delta_row, delta_col]);

FT_centrada6 = FT(Gminus).*circ1;
Centered_complex_field6 = IFT(FT_centrada6);
zero_phase6 = (Centered_complex_field6);
end