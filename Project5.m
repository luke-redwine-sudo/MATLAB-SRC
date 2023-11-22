% Author: Luke Redwine
% Date: 11/20/2022
% Class: ECE 8473 - Digital Image Processing
% Project 5

clc;
clear;
close all;

% Only allow photo files to be selected
file = uigetfile('*.jpg;*.png;*.tif;*.jpeg','Select a File');

% Gray the image
grayed_image = double(rgb2gray(imread(file)));

% Convert the image to a histogram of 1 to 255
raw_hist = hist(grayed_image(:), 255);

% Divide each level in the histogram by the sum of the histogram
% in order to normalize it
for i = 1:255
    normalized_hist(i) = raw_hist(i) / sum(raw_hist);

    % If the histogram value is nan, return 0
    if isnan(normalized_hist(i))
        normalized_hist(i) = 0;
    end
end

% Establish the default Qb max rating
OB_MAX = -1;

% Establish the list of Qbs
OB_COL = [];

% Establish the default Threshold
T = -1;

% For each k, 1 to 1 - the length of the histogram
% Determine P1 and P2 and calculate Qb
for k = 1:254

    P1 = 0;
    P2 = 0;

    % Calculate P1 and P2
    for i = 1:255
        if i < k
            P1 = P1 + normalized_hist(i);
        end

        if i > k
            P2 = P2 + normalized_hist(i);
        end
    end

    M1 = 0;
    M2 = 0;

    % Calculate M1 and M2 from the P1 and P2 values
    for i = 1:255

        % If less than current k, calculate M1
        if i < k
            M1 = M1 + ((i * normalized_hist(i)) / P1);
        end

        % If greater than current k, calculate M2
        if i > k
            M2 = M2 + ((i * normalized_hist(i)) / P2);
        end
    end

    % Calculate Qb
    OB = P1 * P2 * (M1 - M2)^2;

    % Append the Qb value to the other Qb values
    OB_COL = [OB_COL, OB];

    % If the newwest Qb value is the maximum, replace QB_MAX with the new
    % values and set the threshold to the current k
    if OB > OB_MAX
        OB_MAX = OB;
        T = k;
    end

end

% For each pixel in the image, determine if the intensity value is greater
% than or less than the calculated threshold value and assign that pixel to
% be white or black accordingly
for i = 1:size(grayed_image, 1)
    for j = 1:size(grayed_image, 2)
        if grayed_image(i, j) > T
            segmented_image(i, j) = 255;
        else
            segmented_image(i, j) = 0;
        end
    end
end

% Display the histogram of the original image with the calculated threshold
% added on
figure()
stem(1:255, raw_hist);
line([T, T], ylim, 'LineWidth', 1, 'Color', 'r');
xlabel("Gray Level");
ylabel("Frequency")
title("Original Image Histogram");

% Display the original image and the segemented image
figure()
subplot(1,2,1);
imshow(grayed_image, []);
title("Original Gray Scale Image (a)");
subplot(1,2,2);
imshow(segmented_image, []);
title("Segmented Image (b)");

% Display the list of calculated Qbs with the calculated threshold added on
figure()
plot(1:254, OB_COL);
line([T, T], ylim, 'LineWidth', 1, 'Color', 'r');
xlabel("k");
ylabel("Sigma");
title("Sigma Qb vs K");
