%% ECE 8473 Digital Image Processing: Project 2, Part 1
% By Rebecca Garcia
% Dated: 9/20/2022
% Simulate motion blur. Show the original image, degradation kernel image,
% and blurred image
clc; clear all
%% Loading File
pic0 = imread('cat.jpg'); % Read image
picgrey = rgb2gray(pic0); % Coverting to grayscale
pic = double(picgrey);    % Converting from uint8 to double for processing
%% Noise Parameters
m = (length(pic(:,1)));
n = (length(pic(1,:)));
N = m*n;
noise = (255-.01).*rand(m,n) + .0001;% Random gaussian Noise distribution from 0 to 255  
zhat = mean(mean(noise));  % Sample Mean
% Intializing
sigma_eq = zeros(N,1); p = zeros(N,1); P = zeros(m,n);
k = 1;
sigma_eq(k) = (noise(1,1)-zhat(1))^2; % Sample Variance
% PDF of Gaussian Noise Loop
for i = 2:m
    for j = 2:n
        k = k+1;
        sigma_eq(k,1) = (sigma_eq(k-1,1)+(noise(i,j)-zhat)^2); % Sample Variance
    end
end

sigma_sq = (1/(N-1)) * sigma_eq;
        
 for i = 2:m
    for j = 2:n       
        cnst = 1/(sqrt(2*pi*sigma_sq(k,1)));
        z_diff = (noise(i,j)-zhat)^2;
        p(k,1) = cnst * exp( -1*(z_diff)/ (2*sigma_sq(k,1)) );
        P(i,j) = p(k,1);
    end
end
%% images
P1 = (25500-.01)*P+ .0001; 
noisypic = pic+P1;
figure(2);
subplot(1,3,1)
title('Original')
imshow(picgrey)
subplot(1,3,2)
title('Gaussian Noise Dist')
imhist(uint8(P1))
subplot(1,3,3)
imshow(uint8(noisypic))
title('noisypic')