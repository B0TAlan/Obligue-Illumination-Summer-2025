clc;
clear all;
close all;

[f,p] = uigetfile("*.bmp",MultiSelect="off");

file = strcat(p,f);
coplex = Ophase_rec(file,4.65,4.65,.532,2,0);

cf = FT(coplex);

figure; colormap gray;imagesc(log(abs(cf).^2));axis image;

