% Author: Luke Redwine
% Date: 09/17/2022
% Class: ECE 8473 - Digital Image Processing
% Project 1 - Part 3

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

% Get the width and length of the image and double it for the new size
len_m = size(grayed_image, 1)*2;
len_n = size(grayed_image, 2)*2;

doubled = zeros(len_m, len_n);
doubled(1:len_m/2,1:len_n/2) = grayed_image;

% Get the half distances for the distance calculation
dist_p = round(len_m/2);
dist_q = round(len_n/2);

% Cutoff Frequency for the Filter
cutFreq = 60;

% Create the Kernal
kernal = zeros(len_m, len_n);
for i=1:len_m
    for j=1:len_n
        % Calculate the distance from the center
        squared_dist = (i - dist_p).^2 + (j - dist_q).^2;

        % Factor the cut off frequency into the constant for the kernal
        kernal(i,j) = exp(-squared_dist / (2*cutFreq*cutFreq));
    end
end

% Create the view of the filter in the frequency domain
spectrum = fftshift(fft2(kernal));

% Create the DFT of F(u,v) by taking the Fourier Transform of the image
DFT_Fuv = fftshift(fft2(doubled));

% Create the Product of G(u,v) with array multiplication
product_Guv = DFT_Fuv.*kernal;

% Take the Inverse Discrete Fourier Transform
processed_Gp = abs(ifft2(product_Guv));

% Take the upper left corner of the image
final_image = processed_Gp(1:len_m/2,1:len_n/2);

figure();
subplot(1,2,1);
imshow(uint8(grayed_image));
title("Original Image (a)");

subplot(1,2,2);
imshow(uint8(final_image));
title("Gaussian Filtered Image Freq = " + cutFreq + " (b)");

figure();
subplot(1,2,1);
imshow(log(abs(DFT_Fuv)), []);
title("Spectrum of the Image (a)");

subplot(1,2,2);
imshow((abs(kernal)), []);
title("Spectrum of the Filter (b)");
