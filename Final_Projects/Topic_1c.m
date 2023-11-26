close all

% Load the image dataset
%load Sentinel2-data.mat;

%Sentinel2_data = struct2cell(load("Sentinel2-data.mat"));

numSlices = 5;

chosenImage = 5;

% Load the image dataset
load pan_IKONOS.mat;

pan_IKONOS = struct2cell(load("pan_IKONOS.mat"));

% Combine the multispectral and panchromatic data
ms_data = struct2cell(load("ms_IKONOS.mat"));

%ms_data = ms_data{1, 1}; % Assuming that the multispectral data is stored in a field named 'data'
pan_data = load("pan_IKONOS.mat"); % Assuming that the panchromatic data is stored in a field named 'data'

panchromatic_image = pan_IKONOS{1, 1};

% Resize multispectral images to match panchromatic image size
[m, n] = size(panchromatic_image);
multispectral_images = cellfun(@(x) imresize(x, [m, n]), ms_data, 'UniformOutput', false);

% Combine multispectral images into a single matrix
multispectral_matrix = cell2mat(reshape(multispectral_images, 1, 1, []));

% Perform PCA on the multispectral images
mean_multispectral = mean(multispectral_matrix, 3);
centered_multispectral = multispectral_matrix - repmat(mean_multispectral, 1, 1, size(multispectral_matrix, 3));

% Reshape the centered_multispectral to a 2D matrix for PCA
reshaped_multispectral = reshape(centered_multispectral, [], size(centered_multispectral, 3));

%panchromatic_image = upsample_ms(panchromatic_image);

reshaped_multispectral(:, 1) = reshape(panchromatic_image, m*n, 1);

% Compute the covariance matrix using X^T*X
covariance_matrix_multispectral = (reshaped_multispectral' * reshaped_multispectral) / (size(reshaped_multispectral, 1) - 1);

% Perform PCA on the covariance matrix
[eigenvectors_multispectral, eigenvalues_multispectral] = eig(covariance_matrix_multispectral);

% Sort eigenvectors based on eigenvalues in descending order
[eigenvalues_multispectral, sorted_indices_multispectral] = sort(diag(eigenvalues_multispectral), 'descend');
eigenvectors_multispectral = eigenvectors_multispectral(:, sorted_indices_multispectral);

% Number of principal components to keep
num_components = 4;

% Use the first 'num_components' principal components for pansharpening
top_pcs_multispectral = eigenvectors_multispectral(:, 1:num_components);

% Pansharpening using the top principal components
pansharpened_components = reshape(reshaped_multispectral * top_pcs_multispectral, m, n, num_components);
pansharpened_image = sum(pansharpened_components, 4) + panchromatic_image;

NIR = rescale(multispectral_images{1, 1}(:,:,4));
R = rescale(multispectral_images{1, 1}(:,:,3));
G = rescale(multispectral_images{1, 1}(:,:,2));
B = rescale(multispectral_images{1, 1}(:,:,1));

NDVI = normalize(((NIR) - (R)) ./ ((NIR) + (R)));

NIR_P = normalize(pansharpened_components(:,:,4));
R_P = normalize(pansharpened_components(:,:,3));
G_P = normalize(pansharpened_components(:,:,2));
B_P = normalize(panchromatic_image);

NDVI_P = (NIR_P - R_P) ./ (NIR_P + R_P);

figure;
subplot(1, 3, 1);
imshow(cat(3, R, G, B));
title(['Multispectral Color Image']);
subplot(1, 3, 2);
imshow(panchromatic_image, []);
title(['Panchromatic Image']);
subplot(1, 3, 3);
imshow(cat(3, B_P, G_P, R_P), []);
title(['PCA Pansharpened Image']);

figure;
subplot(1, 2, 1);
imshow(cat(3, NDVI, zeros(m, n), zeros(m,n)), []);
title(['NDVI Before Pansharpening'])
subplot(1, 2, 2);
imshow(cat(3, NDVI_P, zeros(m, n), zeros(m,n)), []);
title(['NDVI After Pansharpening']);


% -------------------------------------------------------------------------

multispectral_bands = multispectral_images{1,1};
denominator = multispectral_bands(:,:,1) + multispectral_bands(:,:,2) + multispectral_bands(:,:,3);

fused_band1 = (multispectral_bands(:,:,1).*panchromatic_image)./denominator;
fused_band2 = (multispectral_bands(:,:,2).*panchromatic_image)./denominator;
fused_band3 = (multispectral_bands(:,:,3).*panchromatic_image)./denominator;

fused_image = cat(3, fused_band1, fused_band2, fused_band3)/255;
figure;
imshow(fused_image);
title(['Brovey Method Pansharpened Image'])
