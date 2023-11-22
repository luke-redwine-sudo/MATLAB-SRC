% Author: Luke Redwine
% Date: 12/02/2022
% Class: ECE 8473 - Digital Image Processing
% Bonus Project

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

% Establish the default Threshold
T1 = -2;
T2 = -1;

% Establish the default Qb max rating
OB_MAX = -1;

% Establish the list of Qbs
OB_COL = [];

% For each k1 and k2, 1 to 1 - the length of the histogram
% Determine P1 and P2 and P3 and calculate Qb
for k1 = 1:254
    for k2 = 1:254

        P1 = 0;
        P2 = 0;
        P3 = 0;
    
        % Calculate P1 and P2
        for i = 1:255
            if i < k1
                P1 = P1 + normalized_hist(i);
            end
    
            if i > k1 & i < k2
                P2 = P2 + normalized_hist(i);
            end

            if i > k2
                P3 = P3 + normalized_hist(i);
            end
        end
    
        M1 = 0;
        M2 = 0;
        M3 = 0;
    
        % Calculate M1, M2, M3 from the P1, P2, P3 values
        for i = 1:255
    
            % If less than current k1, calculate M1
            if i < k1
                M1 = M1 + ((i * normalized_hist(i)) / P1);
            end
    
            % If greater than current k1 and less than current k2, calculate M2
            if i > k1 & i < k2
                M2 = M2 + ((i * normalized_hist(i)) / P2);
            end

            % If greater than current k2, calculate M2
            if i > k2
                M3 = M3 + ((i * normalized_hist(i)) / P3);
            end
        end
    
        % Calculate MT
        MT = (P1 * M1) + (P2 * M2) + (P3 * M3);

        % Calculate Qb    
        OB = (P1 * (M1 - MT)^2) + (P2 * (M2 - MT)^2) + (P3 * (M3 - MT)^2);

        % Append the Qb value to the other Qb values
        OB_COL(k1, k2) = OB;
        
        % If the newwest Qb value is the maximum, replace QB_MAX with the new
        % values and set the threshold to the current k1 and k2
        if OB > OB_MAX & k1 < k2
            OB_MAX = OB;
            T1 = k1;
            T2 = k2;
        end
    end
end

% Set all values where K1 is greater than K2 to 0
for i = 1:254
    for j = 1:254
        if i > j
            OB_COL(i, j) = 0;
        end
    end
end

% For each pixel in the image, determine if the intensity value is greater
% than or less than the calculated threshold value and assign that pixel to
% be white or black accordingly
for i = 1:size(grayed_image, 1)
    for j = 1:size(grayed_image, 2)
        if grayed_image(i, j) < T2 & grayed_image(i, j) > T1
            segmented_image(i, j) = 128;
        elseif grayed_image(i, j) > T2
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
line([T1, T1], ylim, 'LineWidth', 1, 'Color', 'r');
text(T1+4,max(raw_hist), "Threshold k1");
line([T2, T2], ylim, 'LineWidth', 1, 'Color', 'r');
text(T2+4,max(raw_hist), "Threshold k2");
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
surf(1:254, 1:254, OB_COL, 'edgecolor', 'none');
hold on
xlabel("k2");
ylabel("k1");
title("Sigma Qb vs K1 vs K2");
plot3(T2,T1,max(max(OB_COL)), '.', color='red', MarkerSize=30)
hold off