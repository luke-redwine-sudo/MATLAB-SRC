clc;
clear;
close all;
% Read the image 
raw_image=imread('photo15697221778714ed8080b7151400x400.jpg');

% Convert to grayscale incase it is color
grayed_image = rgb2gray(raw_image);

figure;
subplot(1,2,1);


grayed_image = im2double(grayed_image);                                   %convert the range of colors from 0-255 to 0-1
[m,n] =size(grayed_image);                                         %Get size of image
      
imhist(grayed_image);
A = fft2(grayed_image);                                        %fourier transform of image
A1 = fftshift(A);                                   %shifting origin
Aabs = abs(A1);                                      %Magnitude of A1 (Frequency domain representation of image)

disp(['The size of the image is ',num2str(m),' x ',num2str(n)]);
d0 = 80;

% Gaussian Low pass filtering

for i=1:m
    for j=1:n
        dgauss(i,j) = (exp(-(i-m/2)^2/(2*d0^2))*exp(-(j-n/2)^2/(2*d0^2)));
    end
end

Bgaussl = dgauss.*A1;
Bgausslmag = dgauss.*Aabs;
Bgaussl1 = fftshift(Bgaussl);
bgaussl = ifft2(Bgaussl1);

subplot(1,2,2);
imhist(bgaussl);

figure;
subplot(1,2,1);
imshow(grayed_image);
title('original image');

subplot(1,2,2);
imshow(uint8(Aabs));
title('Frequency domain image');

figure;
subplot(1,2,1);
imshow(dgauss);
title('Gaussian Low pass Filter');

subplot(1,2,2);
imshow(abs(bgaussl));
title('Gaussian Low pass Filtered image');
