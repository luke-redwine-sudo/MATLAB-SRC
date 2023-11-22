% Author: Luke Redwine
% Date: 09/05/2022
% Class: ECE 8473 - Digital Image Processing
% Project 1 - Part 2

clc
clear all
close all

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
    grayed_image = double(grayed_image);
end

% Establish figure for the images
figure;
subplot(1,2,1);
imshow(grayed_image, []);
title('Gray Image (a)');

% Settings for the Gaussian Kernel
window_size = 3;
sigma = window_size/6;

floor_size = floor(window_size / 2);

% Create the kernal matrix
[x,y] = meshgrid(-floor_size:floor_size, -floor_size:floor_size);

% Take the M and N sizes of the kernal matrix
M = size(x, 1) - 1;
N = size(y, 1) - 1;

% Implement the formula for a Guassian Kernal
formula = exp(-(x.^2 + y.^2) / (2 * sigma * sigma));
formula_prefix = 1 / sum(formula(:));
kernel = formula_prefix * formula;

% Initialize an empty matrix the size of the input image
resultant_image = zeros(size(grayed_image));

% Pad the matrix with zeros the size of the window
padded_image = padarray(grayed_image, [floor_size floor_size]);

% Perform the convolution over the image with the window
for i = 1:size(padded_image,1) - M
    for j =1:size(padded_image, 2) - N
        kernal_window = padded_image(i:i + M,j:j + N).*kernel;
        resultant_image(i,j) = sum(kernal_window(:));
    end
end

% Image with applied Guassian Smoothing
subplot(1,2,2);
imshow(uint8(resultant_image));
title("Gaussian Low Pass Filter ws = " + window_size + " (b)");

% Establish figure for the histograms
% Histogram of the grayed image
figure;
subplot(1,2,1);
imhist(uint8(grayed_image));
title('Gray Image Histogram (a)');
xlabel("Intensity Values");
ylabel("Number of Pixels");

% Histogram of the filtered image
subplot(1,2,2);
imhist(uint8(resultant_image));
title("Gaussian Low Pass Filter Histogram ws = " + window_size + " (b)");
xlabel("Intensity Values");
ylabel("Number of Pixels");
