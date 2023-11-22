% Author: Luke Redwine
% Date: 10/14/2022
% Class: ECE 8473 - Digital Image Processing
% Project 3

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

% Define frequencies
red_freq = 2;
green_freq = 4;
blue_freq = 3;

% Define thetas
red_theta = 5;
green_theta = 10;
blue_theta = 15;

len_m = size(grayed_image, 1);
len_n = size(grayed_image, 2);

red_image = zeros(len_m, len_n);
green_image = zeros(len_m, len_n);
blue_image = zeros(len_m, len_n);

normalized_image = zeros(len_m, len_n);
grayed_image_min = min(grayed_image(:));
grayed_image_max = max(grayed_image(:));

% Normalize the gray scale image
for i = 1:len_m
    for j = 1:len_n
        normalized_image(i, j) = (grayed_image(i, j)-grayed_image_min)/(grayed_image_max - grayed_image_min);
    end
end

% Fill the red, green, and blue image channels with the sinusoidal
% transformation function
for i = 1:len_m
    for j = 1:len_n
        red_image(i, j) = sin(red_freq * normalized_image(i, j) + red_theta);
        green_image(i, j) = sin(green_freq * normalized_image(i, j) + green_theta);
        blue_image(i, j) = sin(blue_freq * normalized_image(i, j) + blue_theta);
    end
end

% Place all 3 channels into their appropriate position
color_image = zeros(len_m, len_n);
color_image(:,:,1) = red_image;
color_image(:,:,2) = green_image;
color_image(:,:,3) = blue_image;

figure 
subplot(1,3,1);
imshow(raw_image);
title("Original Image (a)");

subplot(1,3,2);
imshow(grayed_image, []);
title("Greyed Image (b)")

subplot(1,3,3);
imshow(color_image, []);
title("Pseudo Colored Image (c)")
