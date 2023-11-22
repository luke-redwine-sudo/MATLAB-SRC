%% ECE 8473 Digital Image Processing: Project 2, Part 1
% By Rebecca Garcia
% Dated: 9/22/2022
% Simulate motion blur. Show the original image, degradation kernel image,
% and blurred image
clc; clear all; close all
%% Loading File
pic0 = imread('cat.jpg'); % Read image
picgrey = rgb2gray(pic0); % Coverting to grayscale
pic = double(picgrey);    % Converting from uint8 to double for processing
%% Making Noise
m = (length(pic(:,1)));
n = (length(pic(1,:)));
N = m*n;
noise = (255+255).*rand(m,n) + .0001;   % Random gaussian Noise, 0 to 255  
noisypic = pic+noise;

%% Images of noise on image
figure(1);
subplot(1,3,1)
title('Original')
imshow(picgrey)
subplot(1,3,2)
title('random noise')
imshow(uint8(noise))
subplot(1,3,3)
title('noisypic')
imshow(uint8(noisypic))

%% Motion Blur
x0 = 0.1;
y0 = 0.1;
T = 1; 
t = T;  % setting to T so that a = x0, b = y0
a = x0*(T/t);
b = y0*(T/t);
j = sqrt(-1);
H = zeros(m,n);
% Building Kernel from center
for u = -m/2+1:m/2
    for v = -n/2+1:n/2
         sum = u*a+v*b;
         if sum == 0
            H(u+m/2,v+n/2) = 1;
         else
            H(u+m/2,v+n/2) = ((T*sin(pi*sum))/(pi*sum))*exp(-j*pi*sum);
         end
    end
end

% Applying Motion Blur
pic = fftshift(fft2(pic));
blurrypic = pic.*H;
blurrypic = ((ifft2(ifftshift(blurrypic))));

%% Motion Blur on Image
figure(2);
subplot(1,3,1)
title('Original')
imshow(picgrey)
subplot(1,3,2)
imshow(uint8((abs(H))),[])
title('Degradation Function')
subplot(1,3,3)
imshow(uint8(blurrypic))
title('Blurry Image')

%% Combining Noise and Blur
noisyblur = pic.*H + fftshift(fft2(noise));
figure(3)
imshow(uint8(ifft2(ifftshift(noisyblur))))
title('Noisy and Blurry Image')

%% Part 2, Wiener Filtering
K = .1;
L = zeros(m,n);
G = noisyblur;
H_conj = conj(H);
H_abssq =H_conj.*H;
for u=1:m
    for v=1:n
        L(u,v) = ((H_conj(u,v))/(H_abssq(u,v)+K));
    end
end
F = L.*G;
Fspatial = abs((ifft2((ifftshift((F))))));
%% Restored Image
figure(4)
imshow(uint8(Fspatial))
title("K = " + K)










