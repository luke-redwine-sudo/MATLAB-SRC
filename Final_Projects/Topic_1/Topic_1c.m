close all

% Load the image dataset
load pan_IKONOS.mat;

pan_IKONOS = struct2cell(load("pan_IKONOS.mat"));

% Combine the multispectral and panchromatic data
ms_data = struct2cell(load("ms_IKONOS.mat"));

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

NIR_P = rescale(pansharpened_components(:,:,4));
R_P = rescale(pansharpened_components(:,:,3));
G_P = rescale(pansharpened_components(:,:,2));
B_P = rescale(panchromatic_image);

NDVI_P = (NIR_P - R_P) ./ (NIR_P + R_P);


figure;
subplot(2, 4, 1);
imshow(cat(3, R, zeros(m, n), zeros(m, n)));
title(['Red']);
subplot(2, 4, 2);
imshow(cat(3, zeros(m, n), G, zeros(m, n)));
title(['Green']);
subplot(2, 4, 3);
imshow(cat(3, zeros(m, n), zeros(m, n), B));
title(['Blue']);
subplot(2, 4, 4);
imshow(cat(3, NIR, zeros(m, n), zeros(m, n)));
title(['NIR']);
subplot(2, 4, 5);
imshow(cat(3, R_P, zeros(m, n), zeros(m, n)));
title(['Red']);
subplot(2, 4, 6);
imshow(cat(3, zeros(m, n), G_P, zeros(m, n)));
title(['Green']);
subplot(2, 4, 7);
imshow(cat(3, zeros(m, n), zeros(m, n), B_P));
title(['Blue']);
subplot(2, 4, 8);
imshow(cat(3, NIR_P, zeros(m, n), zeros(m, n)));
title(['NIR']);

figure;
subplot(1, 3, 1);
imshow(uint8(rescale(cat(3, R, G, B), 0, 255)));
title(['Multispectral Color Image']);
subplot(1, 3, 2);
imshow(panchromatic_image, []);
title(['Panchromatic Image']);
subplot(1, 3, 3);
imshow((uint8(rescale(cat(3, R_P, B_P, G_P), 0, 255))));
title(['PCA Pansharpened Image']);

figure;
subplot(1, 2, 1);
imshow(cat(3, NDVI, zeros(m, n), zeros(m,n)), []);
title(['NDVI Before Pansharpening'])
subplot(1, 2, 2);
imshow(cat(3, NDVI_P, zeros(m, n), zeros(m,n)));
title(['NDVI After Pansharpening']);


% -------------------------------------------------------------------------

multispectral_bands = multispectral_images{1,1};
denominator = multispectral_bands(:,:,1) + multispectral_bands(:,:,2) + multispectral_bands(:,:,3);

fused_band1 = (multispectral_bands(:,:,1).*panchromatic_image)./denominator;
fused_band2 = (multispectral_bands(:,:,2).*panchromatic_image)./denominator;
fused_band3 = (multispectral_bands(:,:,3).*panchromatic_image)./denominator;

fused_image = uint8(rescale(cat(3, fused_band1, fused_band2, fused_band3), 0, 255));
figure;
imshow(fused_image);
title(['Brovey Method Pansharpened Image'])
