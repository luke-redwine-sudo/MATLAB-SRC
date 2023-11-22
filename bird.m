a = imread('bird.jpg'); % read an image
figure; image(a); % display a color image

% change a color image to a gray-scale image
a = double(a);
a = (a(:,:,1)+a(:,:,2)+a(:,:,3))/3;
figure; imshow(a,[]) % diplay a gray-scale image
imwrite(a, 'bird_gray.jpg'); % save an image
title('Original image')

% Pixel-level transformation
a = (a-min(a(:)))/(max(a(:))-min(a(:)));
b1 = 1./(1+exp(-5*(2*a-1)));
figure; subplot(1,3,1); imshow(b1,[]); title('After Transform 1')

b2 = a.^2;
subplot(1,3,2); imshow(b2,[]); title('After Transform 2')

b3 = log(10*a+1); b3=b3/max(b3(:));
subplot(1,3,3); imshow(b3,[]); title('After Transform 3')

% Transform functions 
x = [0:0.01:1];
y1 = 1./(1+exp(-5*(2*x-1))); y1 = (y1-min(y1))/(max(y1)-min(y1));
figure; subplot(1,3,1); plot(x,y1); title('Transform function 1')
y2 = x.^2;
subplot(1,3,2); plot(x,y2); title('Transform function 2')
y3 = log(10*x+1); y3 = y3/max(y3(:));
subplot(1,3,3); plot(x,y3); title('Transform function 3')
