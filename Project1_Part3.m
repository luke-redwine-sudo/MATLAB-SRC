% Author: Luke Redwine
% Date: 09/08/2022
% Class: ECE 8473 - Digital Image Processing
% Project 1 - Part 3

clc;
clear;
close all;

% Only allow jpg files to be selected
file = uigetfile('*.jpg;*.png','Select a File');

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
    disp("Image must be grayed")
    grayed_image = (0.2989 *raw_image(:,:,1) + 0.5870 *raw_image(:,:,2) + 0.1140 *raw_image(:,:,3));
    grayed_image = double(grayed_image) / 255;
end

% Get the size of image
[m,n] = size(grayed_image);

% Get the fourier transform of image
% Shift the origin of the image
fourtrans = fft2(grayed_image);
fourtransshift = fftshift(fourtrans);

% Obtain the Magnitude of the Shifted Fourier Transform to reach the
% Frequency Domain
fourtransshiftabs = abs(fourtransshift);

% Establish the cutoff frequency
cutfreq = 60;

% Gaussian Low pass filtering
gausslowpass = 0;

% Create the guassian low pass filter
for i=1:m
    for j=1:n
        gausslowpass(i,j) = exp(-(i-m/2)^2/(2*cutfreq^2));
        gausslowpass(i,j) = gausslowpass(i,j) * exp(-(j-n/2)^2/(2*cutfreq^2));
    end
end

% Apply the guassian filter to the fourier transform shifted image
filter_resultant = gausslowpass.*fourtransshift;

% Shift the origin point back to its original location
filter_resultant = fftshift(filter_resultant);

% Invert the fourier transform
filter_resultant = ifft2(filter_resultant);

% Create the figure for the images
figure;
subplot(1,2,1);
imshow(grayed_image);
title('Gray Image (a)');

subplot(1,2,2);
imshow(filter_resultant);
title('Gaussian Low pass Filtered Image (b)');

% Create the figure for the histograms
figure;
subplot(1,2,1);
imhist(grayed_image);
title("Gray Image Histogram (c)");

subplot(1,2,2);
imhist(filter_resultant);
title("Filtered Image Histogram (d)")

% Create the figure for the filters
figure;
subplot(1,2,1);
imshow(gausslowpass);
title('Gaussian Low Pass Filter (e)');

subplot(1,2,2);
imshow(uint8(fourtransshiftabs));
title('Frequency Domain Image (f)');
