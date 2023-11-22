% Author: Luke Redwine
% Date: 09/22/2022
% Class: ECE 8473 - Digital Image Processing
% Project 2 - Part 2

clc;
clear;
close all;

% Only allow files to be selected
file = uigetfile('*.jpg;*.png;*.tif','Select a File');

% Determine if a file was selected
% If it was not, exit the code
if isequal(file,0)
    disp('File required for execution. Exiting...');
    return
else
    disp('File Selected');
    raw_image = imread(file);
end

if size(raw_image, 3) == 1
    disp("Image is already grayed")
    grayed_image = raw_image;
else
    disp("Image must be grayed");
    grayed_image = double(rgb2gray(raw_image));
end


len_m = size(grayed_image, 1);
len_n = size(grayed_image, 2);

G = fftshift(fft2(grayed_image));

H = idft(grayed_image);

K = 0.025;
L = zeros(len_m, len_n);
for u = 1:len_m
    for v = 1:len_n
        L(u, v) = conj(H(u, v))/((conj(H(u, v))*H(u, v)) + K);
    end
end

F = ifft2(fftshift(G.*L));

figure
imshow(log(abs(L)), []);

figure
imshow(real(F), []);