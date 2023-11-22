clc;
clear;
close all;

a = imread('book.tif'); figure; imshow(a,[]);
grayed_image = im2gray(a);

% Get the size of image
[m,n] = size(grayed_image);
p = m * 2;
q = n * 2;

% Create a new matrix array and place the image in the upper left corner
a = zeros(p,q);
a(1:m,1:n) = grayed_image;

%a = double(grayed_image);
ya = fft2(a);
figure; imshow(abs(ya),[])
figure; imshow(log(abs(ya)),[])

a2 = a;
for i = 1:size(a,1)
    for j = 1:size(a,2)
        a2(i,j) = a(i,j)*(-1)^(i+j);
    end
end
ya2 = fft2(a2);
figure; imshow(abs(ya2),[])
figure; imshow(log(abs(ya2)),[])

ya3=fftshift(fft2(a)); figure; imshow(log(abs(ya3)),[])
