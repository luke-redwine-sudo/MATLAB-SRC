% Author: Luke Redwine
% Date: 09/05/2022
% Class: ECE 8473 - Digital Image Processing
% Project 1 - Part 1

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
    grayed_image = double(raw_image);
    grayed_image = (0.2989 *raw_image(:,:,1) + 0.5870 *raw_image(:,:,2) + 0.1140 *raw_image(:,:,3));
end

% Establish figure for the images
figure;

% Plot for pre-equalization image
subplot(2,2,1);
imshow(grayed_image);
title('Gray Image (a)');

% Plot for pre-equalization histogram
subplot(2,2,3);
imhist(grayed_image);
title('Gray Image Histogram (c)');
xlabel("Intensity Values");
ylabel("Number of Pixels");

% Gather the rows and columns of the grayed image
rows = size(grayed_image, 1);
columns = size(grayed_image, 2);

% Get the total number of pixels in the image
pixels = rows * columns;

% Create blank arrays for frequency of occurance, pdf, cdf, and sk function
frequency = zeros(256,1);
pdf = zeros(256,1);
cdf = zeros(256,1);
sk_function = zeros(256,1);

% Examine the intensity of the pixels individually. Use the frequency of
% occurance over the total number of pixels to determine the pdf
for i = 1:rows
    for j = 1:columns
        pixel_value = grayed_image(i,j);
        % Add one to the pixel value so no value is zero
        frequency(pixel_value+1) = frequency(pixel_value+1)+1;
        pdf(pixel_value+1) = frequency(pixel_value+1)/pixels;
    end
end

sum = 0;

% Apply the pdf over the frequency array to create the cdf. Divide the
% total number of pixels traversed over the total number of pixels to
% create the CDF. Then multiply the value at i by the light levels and
% round so it displays correctly
for i = 1:size(pdf)
    sum = sum + frequency(i);
    cdf(i) = sum/pixels;
    sk_function(i) = round(cdf(i)*256);
end

% Apply the acquired sk function to the image
for i = 1:rows
    for j = 1:columns
        resultant_image(i, j) = sk_function(grayed_image(i,j)+1);
    end
end

% Plot for post-equalization image
subplot(2,2,2);
imshow(uint8(resultant_image));
title("Resultant Image (b)");

% Plot for post-equalization histogram
subplot(2,2,4);
imhist(uint8(resultant_image));
title('Resultant Image Histogram (d)');
xlabel("Intensity Values");
ylabel("Number of Pixels");
