% Author: Luke Redwine
% Date: 09/18/2022
% Class: ECE 8473 - Digital Image Processing
% Project 2 - Part 1

clc;
clear;
close all;

% Only allow jpg files to be selected
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
len_m = size(grayed_image, 1);
len_n = size(grayed_image, 2);

num_pixels = len_m * len_n;

% Calculate the z mean value
z_mean = sum(grayed_image, "all") / num_pixels;


% Calculate the variance for the guassian noise filter
variance = 0;
for i=1:len_m
    for j=1:len_n
        pixel = grayed_image(i,j);
        variance = variance + ((grayed_image(i,j) - z_mean)^2);
    end
end

% Implement Variance formula
variance = variance / (num_pixels - 1);

% Create the Guassian noise constant from the book
const_form = 1 / (sqrt(2 * pi * variance));

% Create a noise filter the same size as the image
noise = zeros(len_m, len_n);
for i=1:len_m
    for j=1:len_n
        noise(i,j) = const_form * exp(-((grayed_image(i,j) - z_mean)^2) / (2 * variance));        
    end
end

% Create a Discrete fourier transform of the image
DFT = fftshift(fft2(grayed_image));

% Apply the gaussian filter and return from the frequency domain
noisy_image = ifft2(ifftshift(DFT.*noise));

figure 
subplot(1,3,1);
imshow(uint8(grayed_image));

subplot(1,3,2)
imshow(abs(noise), []);

subplot(1,3,3)
imshow(abs(noisy_image), []);