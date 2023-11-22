% Author: Luke Redwine
% Date: 09/08/2022
% Class: ECE 8473 - Digital Image Processing
% Project 1 - Part 3

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
    grayed_image = (0.2989 *raw_image(:,:,1) + 0.5870 *raw_image(:,:,2) + 0.1140 *raw_image(:,:,3));
    grayed_image = double(grayed_image) / 255;
end

% Settings for the Gaussian Kernel
window_size = 3;
sigma = window_size/6;

floor_size = floor(window_size / 2);

% Create the kernal matrix
[x,y] = meshgrid(-floor_size:floor_size, -floor_size:floor_size);

% Take the M and N sizes of the kernal matrix
M = size(x, 1);
N = size(y, 1);

% Get the size of image
[m,n] = size(grayed_image);
p = m * 2;
q = n * 2;

% Create a new matrix array and place the image in the upper left corner
doubled_image = zeros(p,q);
doubled_image(1:m,1:n) = grayed_image;

% Implement the formula for a Guassian Kernal
formula = exp(-(x.^2 + y.^2) / (2 * sigma * sigma));
formula_prefix = 1 / sum(formula(:));
kernal = formula_prefix * formula;
padded_kernal = zeros(length(doubled_image(:,1)), length(doubled_image(1,:)));
padded_kernal(1:M, 1:N) = kernal;

%padded_kernal = fftshift(fft2(padded_kernal));

% Implement the formula for a Guassian Kernal
% formula = zeros(M,N);
% kernal = formula;
% cutoffFreq = 20;
% for i = 1:M
%     for j = 1:N
%         formula(i,j) = sqrt(((i-(M/2))^2+(j-(N/2))^2));
%         kernal(i,j) = exp(-(formula(i,j)^2)/(2*cutoffFreq^2));
%      end
% end

full_filter = zeros(p,q);

for i = 1:m
    for j = 1:n
        full_filter(i,j) = doubled_image(i,j)*((-1)^(i+j));
    end
end

dft_filter = fftshift(fft2(full_filter));

% % Perform the convolution over the image with the window
% padded_filter = padarray(dft_filter, [floor_size floor_size]);
% 
% H = 0;
% for i = 1:size(padded_filter,1) - M + 1
%     for j =1:size(padded_filter, 2) - N + 1
%         kernal_window = padded_filter(i:i + M - 1,j:j + N - 1).*kernal;
%         H(i,j) = sum(kernal_window(:));
%         %H(i:i + M - 1,j:j + N - 1) = padded_filter(i:i + M - 1,j:j + N - 1).*kernal;
%     end
% end

guassian_product = zeros(length(doubled_image(:,1)), length(doubled_image(1,:)));

for i=1:p
    for j=1:q
        temp = full_filter(i,j) * kernal(i,j);
        guassian_product(i,j) = temp;
    end
end

%guassian_product = full_filter.*kernal;
guassian_product = ifft2(ifftshift(guassian_product));

processed_image = zeros(p,q);

for i=1:length(guassian_product(:,1))
    for j=1:length(guassian_product(1,:))
        processed_image(i,j) = real(guassian_product(i,j))*((-1)^(i+j));
    end
end

g = processed_image(1:m,1:n);

figure;
imshow(grayed_image);
title("Grayed Image");

figure;
imshow(doubled_image, []);
title("Doubled Image");

figure;
imshow(log(abs(dft_filter)),[]);
title("Spectrum of Fp")

figure();
imshow(log(abs(H)),[]);
title("Centered Guassian Filter")

figure();
imshow(log(abs(guassian_product)),[]);
title("G(u,v)")

figure();
imshow(log(abs(processed_image)),[]);
title("Gp(u,v)");

figure();
imshow(uint8(g),[]);
title("Top Left Image");
