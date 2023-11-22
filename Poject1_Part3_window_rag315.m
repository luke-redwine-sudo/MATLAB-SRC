%% PART 3
clc;
clear all;
close all;
% PreProccessing Image
pic0 = imread('Screenshot 2022-09-15 232411.jpg');
pic = rgb2gray(pic0);
pic = double(pic);
M = length(pic(:,1));
N = length(pic(1,:));
% Taking FFt
Pic = fftshift(fft2(pic));
%% Building Window
m = 133; n = 134;

p = round(m/2);
q = round(n/2);

sigma = 10; 
kernel = zeros(m,n);
for i=1:m
  for j=1:n
      sq_dist = (i-p).^2 + (j-q).^2;
      kernel(i,j) = exp(-sq_dist/(2*sigma^2));
  end
end
p = abs(length(kernel(:,1))-M);
q = abs(length(kernel(1,:))-N);
A = padarray(kernel, [p/2 q/2]);
H = fftshift(fft2(A));


% Apply the gaussian filter
B = Pic.*A;
C = abs(ifft2(B));

% now we can show original and filtered image
figure(1);
set(gcf, 'Position', get(0, 'Screensize'));
subplot(121),imshow(pic0),title('Original Image');
subplot(122), imshow(uint8(C)), title('output of gaussian filter');