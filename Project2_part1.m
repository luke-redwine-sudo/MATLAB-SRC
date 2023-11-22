% Author: Luke Redwine
% Date: 09/18/2022
% Class: ECE 8473 - Digital Image Processing
% Project 2

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

% Send the image and the A, B, and T values for blurring
[blurred_image, blurring_filter] = blur_image(grayed_image, 0.1, 0.1, 1);

% Send the blurred image with the desired variance value
[noisy_image, noise] = guassian_noise(blurred_image, 3);

% Send the blurred image 
[repaired_image, wiener_filter] = fixing_filter(noisy_image, blurring_filter, 0.001);

figure 
subplot(1,2,1);
imshow(uint8(grayed_image));
title("Original Image (a)");

subplot(1,2,2)
imshow(abs(blurring_filter), []);
title("Motion Blur Kernal (a)");

figure
subplot(1,2,1)
imshow(real(blurred_image), []);
title("Image with Motion Blur Kernal Applied (a)");

subplot(1,2,2)
imshow(real(noisy_image), []);
title("Image with Motion Blur and Noise (b)")

figure 
subplot(1,2,1);
imshow(uint8(grayed_image));
title("Original Image (a)");

subplot(1,2,2)
imshow(real(noisy_image), []);
title("Image with Motion Blur and Noise (b)")

figure 
subplot(1,2,1);
imshow(abs(wiener_filter), []);
title("Wiener Filter (a)");

subplot(1,2,2)
imshow(real(repaired_image), []);
title("Deblurred Image (b)");

figure 
subplot(1,2,1);
imshow(uint8(grayed_image));
title("Original Image (a)");

subplot(1,2,2);
imshow(real(repaired_image), []);
title("Deblurred Image (b)")

figure 
subplot(1,2,1);
imshow(real(noisy_image), []);
title("Image with Motion Blur and Noise (a)")

subplot(1,2,2);
imshow(real(repaired_image), []);
title("Deblurred Image (b)")

function [blurred_image, blurring_filter] = blur_image(image, a, b, T)
    % Get the width and length of the image and double it for the new size
    len_m = size(image, 1);
    len_n = size(image, 2);

    % Create an emtpy matrix for the blurring kernal
    blurring_filter = zeros(len_m, len_n);

    % Place in the image into the frequency domain
    DFT = fftshift(fft2(image));

    % Create the Motion Blur kernal
    for u=-len_m/2+1:len_m/2
        for v=-len_n/2+1:len_n/2
            if (u*a + v*b) == 0
                blurring_filter(u+len_m/2, v+len_n/2) = 1;
            else
                blurring_filter(u+len_m/2, v+len_n/2) = blurring_formula(u, v, T, a, b, j);
            end
        end
    end

    % Return the image in the spatial domain with the blurred filtered
    % applied for display purposes
    blurred_image = ifft2(fftshift(DFT.*blurring_filter));

end

function [noisy_image, noise] = guassian_noise(image, variance)
    % Get the width and length of the image and double it for the new size
    len_m = size(image, 1);
    len_n = size(image, 2);
    
    % Create the noise kernal
    noise = variance * randn(len_m, len_n);
    noise = fftshift(fft2(noise));
    
    % Add the noise kernal to the frequency domain image and convert back
    % to the spatial domain
    noisy_image = ifft2(fftshift(fftshift(fft2(image)) + noise));   
    
end

function blur_output = blurring_formula(i, j, T, a, b, complex_number)
    % Create the different parts of the sinusoidal function
    prefix = (T/(pi*(i*a + j*b)));
    formula = sin(pi*(i*a + j*b));
    suffix = exp(-(complex_number*pi*(i*a + j*b)));

    % Combine the parts of the formula
    blur_output = prefix * formula * suffix;
end

function [repaired_image, wiener_filter] = fixing_filter(blurred_image, degradation_function, K)
    % Get the width and length of the image and double it for the new size
    len_m = size(blurred_image, 1);
    len_n = size(blurred_image, 2);
    
    % Get the fourier transform of the degraded image
    G = fftshift(fft2(blurred_image));

    % Get the complex part of the degradation transfer function
    degradation_conj = conj(degradation_function);
    degradation_abssq = degradation_function.*degradation_conj;

    % Create the Wiener Filter
    for u=1:len_m
        for v=1:len_n
            wiener_filter(u, v) = degradation_conj(u,v)/(degradation_abssq(u,v)+K);
        end
    end
    
    % Apply the Wiener Filter to the blurred image in the frequency domain
    repaired_image = wiener_filter.*G;
    
    % Transfer the image back to the spatial domain
    repaired_image = abs(ifft2(ifftshift(repaired_image)));
end