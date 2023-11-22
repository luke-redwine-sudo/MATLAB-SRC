%% restart
close all;
clc;

%% read Image and Display
test = imread('photo15697221778714ed8080b7151400x400.jpg'); %read image
figure,imshow(test);
figure, histogram(test);

test = double(test);
test = uint8((test(:,:,1) + test(:,:,2) + test(:,:,3)) / 3);


%% Histogram equalization without matlab function
I = test; 

% size of the input image.
[r,c] = size(I); 

% uint 8 blamk canvas
blank = uint8(zeros(r,c));

% number of pixels.
n = r*c;

% frequency pdf and cdf and some variables
f = zeros(256,1);
pdf = zeros(256,1);
cdf = zeros(256,1);
out = zeros(256,1);
cum = zeros(256,1);


% Nested for loop

for i = 1:r
    for j = 1:c
        value = I(i,j);
        f(value+1) = f(value+1)+1;
        pdf(value+1) = f(value+1)/n;
    end
end

% finding cdf
sum = 0;
L = 255;

for i = 1:size(pdf);
    sum = sum + f(i);
    cum(i) = sum;
    cdf(i) = cum(i)/n;
    out(i) = round(cdf(i)*L);
end

for i = 1:r;
    for j = 1:c;
        blank(i,j) = out(I(i,j)+1);
    end
end

figure,imshow(blank); title('My hist Image');
figure, histogram(blank); title('My funtion Histogram');