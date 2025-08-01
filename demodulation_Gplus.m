function [circ1,zero_phase6] = demodulation_Gplus(std_temp_gplus, angles,FTHolo, dxy, lambda, circ2)
posicion = find(std_temp_gplus == min(std_temp_gplus));
theta_final = [0, angles(posicion)];
[Gdemod] =  demComp2SIDHM(theta_final, FTHolo);
Gplus_demod = Gdemod(:,:,1);
gplus = fftshift(ifft2(fftshift(Gplus_demod)));
gplus = gplus*-1;
[ref_gplus] = phase_rec(gplus, dxy, dxy, lambda, circ2);
Gplus = gplus .* ref_gplus;

Gminus1 = FT(Gplus);
% Find max in A_before
[~, idx_before] = max(log(abs(Gplus_demod(:)).^2));
[row_before, col_before] = ind2sub(size(Gplus_demod), idx_before);

% Find max in A_after
[~, idx_after] = (max(log(abs(Gminus1(:)).^2)));
[row_after, col_after] = ind2sub(size(Gminus1), idx_after);

% Calculate displacement
delta_row = row_after - row_before;
delta_col = col_after - col_before;
circ1 = circshift(circ2, [delta_row, delta_col]);

FT_centrada6 = FT(Gplus).*circ1;
Centered_complex_field6 = IFT(FT_centrada6);
zero_phase6 = (Centered_complex_field6);
end