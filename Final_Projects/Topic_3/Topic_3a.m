close all
clear all

load("Sentinel2-data.mat");

%sentinel_6 = {B2, B3, B4, B8, B11, B12};
sentinel_6 = {B2, B3, B4};
sentinel_13 = {B1, B2, B3, B4, B5, B6, B7, B8, B8A, B9, B10, B11, B12};

resize_rows = 1562;
resize_columns = 1524;
resize_factor = [resize_rows resize_columns];

for B = 1:length(sentinel_6)
    sentinel_6{B} = imresize(sentinel_6{B}, resize_factor);
end

for B = 1:length(sentinel_13)
    sentinel_13{B} = imresize(sentinel_13{B}, resize_factor);
end

% Define  Tasseled Cap coefficients
% sentinel_coeffs = [0.3510, 0.3813, 0.3437, 0.7196, 0.2396, 0.1949;
%                   -0.3599, -0.3533, -0.4734, 0.6633, 0.0087, -0.2856;
%                   0.2578, 0.2305, 0.0883, 0.1071, -0.7611, -0.5308;
%                   0.0805, - 0.0498, 0.1950, -0.1327, 0.5752, -0.7775;
%                   -0.7252, -0.0202, 0.6683, 0.0631, -0.1494, -0.0274;
%                   0.400, -0.8172, 0.3832, 0.0602, -0.1095, 0.0985];

sentinel_coeffs = [0.141, 0.019, 0.325;
                   -0.016, 0.629, 0.154;
                   0.714, 0.315, -0.03]

bands = 3;
reshaped_image = cat(3, sentinel_6{:}); % Concatenate along the third dimension
reshaped_image = reshape(reshaped_image, resize_rows * resize_columns, bands);

% Perform Tasseled Cap Transformation
tct_image = double(reshaped_image) * sentinel_coeffs';

% Reshape the result back to the original image dimensions
tct_image = reshape(tct_image, resize_rows, resize_columns, size(sentinel_coeffs, 1));

%SWIR1 = normalize(tct_image(:,:,6));
%SWIR2 = normalize(tct_image(:,:,5));
%NIR = normalize(tct_image(:,:,4));
Red = normalize(tct_image(:,:,3));
Green = normalize(tct_image(:,:,2));
Blue = normalize(tct_image(:,:,1));

%NDVI = ((NIR - Red) ./ (NIR + Red));
%NDWI = ((Green - NIR) ./ (Green + NIR));
%NDSI = ((Green - SWIR1) ./ (Green + SWIR1));

figure;
color = uint8(rescale(cat(3, Red, Green, Blue), 0, 255));
imshow(color);
title("Color");

figure;
subplot(1,3,1);
imshow(cat(3, Red, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("Red");
subplot(1,3,2);
imshow(cat(3, NIR, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("NIR");
subplot(1,3,3);
imshow(cat(3, NDVI, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("6 Band NDVI");

figure;
subplot(1,3,1);
imshow(cat(3, zeros(resize_rows, resize_columns), Green, zeros(resize_rows, resize_columns)), []);
title("Green");
subplot(1,3,2);
imshow(cat(3, NIR, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("NIR");
subplot(1,3,3);
imshow(cat(3, NDWI, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("6 Band NDWI");

figure;
subplot(1,3,1);
imshow(cat(3, zeros(resize_rows, resize_columns), Green, zeros(resize_rows, resize_columns)), []);
title("Green");
subplot(1,3,2);
imshow(cat(3, SWIR1, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("SWIR1");
subplot(1,3,3);
imshow(cat(3, NDSI, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("6 Band NDSI");

figure;
subplot(1,3,1);
imshow(cat(3, NDVI, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("6 Band NDVI");
subplot(1,3,2);
imshow(cat(3, NDWI, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("6 Band NDWI");
subplot(1,3,3);
imshow(cat(3, NDSI, zeros(resize_rows, resize_columns), zeros(resize_rows, resize_columns)), []);
title("6 Band NDSI");