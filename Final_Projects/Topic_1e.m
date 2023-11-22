% Author: Luke Redwine
% Date: 10/16/2022
% Class: ECE 8473 - Digital Image Processing
% Final Project

clc;
clear;
close all;

ms_data = struct2cell(load("ms_IKONOS.mat", "-mat"));
ms_data = ms_data{1, 1};

pan_data = struct2cell(load("pan_IKONOS.mat", "-mat"));
pan_data = pan_data{1, 1};

% Extrract the different dimensions of the struct
lenL = size(ms_data, 3);
lenM = size(ms_data, 1);
lenN = size(ms_data, 2);

% Reshape the data to a 2D array
ms_data = reshape(ms_data, [lenM*lenN lenL])';

% Collect the mean of the ms and then transpose it
meanMs = mean(ms_data')';

% Collect the covariance matrix
covarianceMatrix = (ms_data * ms_data') / (lenM * lenN) - meanMs * meanMs';

% Retrieve the eigenvector values
[eigenVector, eigenValue] = eig(covarianceMatrix);

% Create the variables to store the top eigen values
firstEigen = [0,0];
secondEigen = [0,0];
thirdEigen = [0,0];

% Find the top eigen values
for i = 1:size(eigenValue, 1)
    evalEigen = eigenValue(i, i);

    if evalEigen > firstEigen(1)
        thirdEigen(1) = secondEigen(1);
        thirdEigen(2) = secondEigen(2);
        secondEigen(1) = firstEigen(1);
        secondEigen(2) = firstEigen(2);
        firstEigen(1) = evalEigen;
        firstEigen(2) = i;
    elseif evalEigen > secondEigen(1)
        thirdEigen(1) = secondEigen(1);
        thirdEigen(2) = secondEigen(2);
        secondEigen(1) = evalEigen;
        secondEigen(2) = i;
    elseif evalEigen > thirdEigen(1)
        thirdEigen(1) = evalEigen;
        thirdEigen(2) = i;
    end
end

% Collect the top eigen vectors to complete pca
v1 = eigenVector(:,firstEigen(2));
v2 = eigenVector(:,secondEigen(2));
v3 = eigenVector(:,thirdEigen(2));

% Multiply the transpose of the eigen vectors by the difference of the
% original ms and the meanMs and reshape it to size M x N
y1 = reshape(v1' * (ms_data - meanMs), [lenM lenN]);
y2 = reshape(v2' * (ms_data - meanMs), [lenM lenN]);
y3 = reshape(v3' * (ms_data - meanMs), [lenM lenN]);

% Plot the independant RGB values
figure
subplot(1,3,1);
imshow(y3, []);
title("R Channel")
subplot(1,3,2);
imshow(y2, []);
title("G Channel")
subplot(1,3,3);
imshow(pan_data, []);
title("B Channel")

% Construct the RGB image with the top PCAs
combinedImage = zeros(lenM, lenN);
combinedImage(:,:,1) = pan_data;
combinedImage(:,:,2) = y2;
combinedImage(:,:,3) = y3;

% Show the RGB
figure
imshow(combinedImage, [])
title("Pseudo Color Image")
