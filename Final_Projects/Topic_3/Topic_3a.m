close all

load("Sentinel2-data.mat");

sentinel_6 = {B1, B2, B3, B4, B5, B6};
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
landsat_coeffs = [0.3510, 0.3813, 0.3437, 0.7196, 0.2396, 0.1949;
                  -0.3599, -0.3533, -0.4734, 0.6633, 0.0087, -0.2856;
                  0.2578, 0.2305, 0.0883, 0.1071, -0.7611, -0.5308];

bands = 6;
reshaped_image = reshape(sentinel_6, resize_rows * resize_columns, bands);

% Perform Tasseled Cap Transformation
tct_image = reshaped_image * landsat_coeffs';

% Reshape the result back to the original image dimensions
tct_image = reshape(tct_image, resize_rows, resize_columns, size(landsat_coeffs, 1));

% Display the TCT result for each band
for i = 1:size(landsat_coeffs, 1)
    figure;
    imshow(tct_image(:,:,i), []);
    title(['TCT Band ', num2str(i)]);
end