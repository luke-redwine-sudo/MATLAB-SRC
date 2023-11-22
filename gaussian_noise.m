
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

blurring_filter = zeros(len_m, len_n);

T=1;
a=0.1;
b=0.1;
variance=2;

DFT_image = fftshift(fft2(grayed_image));
noise = fftshift(fft2(variance * randn(len_m, len_n)));

for u=-len_m/2+1:len_m/2
    for v=-len_n/2+1:len_n/2
        prefix = (T/(pi*(u*a + v*b)));
        formula = sin(pi*(u*a + v*b));
        suffix = exp(-(j*pi*(u*a + v*b)));

        if (u*a + v*b) == 0
            blurring_filter(u+len_m/2, v+len_n/2) = 1;
        else
            blurring_filter(u+len_m/2, v+len_n/2) = prefix * formula * suffix;
        end

        
    end
end

blurred_image = DFT_image.*blurring_filter;

final_image = ifft2(fftshift(blurred_image + noise));

figure;
imshow(log(abs(blurred_image)), []);

figure;
imshow(real(final_image), []);