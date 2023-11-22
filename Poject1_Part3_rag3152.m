%% Part 3
% Implement Gaussian Lowpass Filtering in the frequency domain
% Be able to change the window size
% Plot the image and filter in the frequency domain
% following algorithm found on pg 312 of the book

clc;
clear;
close all;

pic0 = imread('Screenshot 2022-09-15 232411.jpg');       % Read image
picgrey = rgb2gray(pic0);        % Coverting to grayscale
pic = double(picgrey);           % Converting from uint8 to double for processing

M = length(pic(:,1));
N = length(pic(1,:));
P = M*2;
Q = N*2;
pic = double(picgrey);

PicPadP3 = zeros(P,Q);
PicPadP3(1:M,1:N) = pic;


figure;
imshow(PicPadP3,[])
PicCenter = zeros(length(PicPadP3(:,1)),length(PicPadP3(1,:)));


for i = 1:P
    for j = 1:Q
        PicCenter(i,j) = PicPadP3(i,j)*(-1)^(i+j);
    end
end
figure;
imshow(PicCenter,[])
title('-1*(x+y)');

PicCenter = fft2(PicCenter);
figure;
imshow(PicCenter,[]);
% get axes limits in pixels
set(gca,'units','pixels');
pos = get(gca,'position');
% display the top left part of the image at magnification 100%
xlim([0.5 pos(3)-0.5]),ylim([0.5 pos(4)-0.5]);
title('fft2');

figure; imshow(log(abs(PicCenter)),[])
title('fft2, log');

% filter window
D = zeros(P,Q);H = D;
D0 = 30;
for i = 1:P
    for j = 1:Q
        D(i,j) = sqrt(((i-(P/2))^2+(j-(Q/2))^2)); %Distance from center
        H(i,j) = exp( -(D(i,j)^2) / (2*D0^2) ); %Gaussian Filter
    end
end
figure;
imshow(H,[]);
title('Filter Window');

G = H.*PicCenter;
G = ifft2(G);
for i=1:length(G(:,1))
    for j=1:length(G(1,:))
            G(i,j) = (real(G(i,j)))*(-1)^(i+j);
    end
end
figure;
imshow(G,[]);
title('g_p');
g = G(1:M,1:N);
figure;
imshow(g,[]);
title('g');
