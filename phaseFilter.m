function [filt, circ] = phaseFilter(file, dx, dy, lambda, region)
im = imread(file);
figure(8), imagesc(im), colormap gray, title 'Image preview'

complexFeild  = Ophase_rec(file, dx,dy,lambda,region);

zero_Phase = phaseNorm(complexFeild);
thickness(zero_Phase, lambda,1.52,1,1);
[f, c] = circular_filter(complexFeild,zero_Phase,12);
filt = f;
circ = c;
end