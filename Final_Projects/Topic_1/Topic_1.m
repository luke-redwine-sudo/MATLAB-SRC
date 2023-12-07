close all
clear all

% Load the image dataset
load pan_IKONOS.mat;

pan_IKONOS = struct2cell(load("pan_IKONOS.mat"));
ms_data = struct2cell(load("ms_IKONOS.mat"));

panchromatic_image = pan_IKONOS{1, 1};

% Retieve the size statistics of the images and retrieve the multspectral
% data
[m, n] = size(panchromatic_image);
multispectral_images = cellfun(@(x) imresize(x, [m, n]), ms_data, 'UniformOutput', false);

% Combine multispectral images into a single matrix
multispectral_matrix = cell2mat(reshape(multispectral_images, 1, 1, []));

% Perform PCA on the multispectral images
mean_multispectral = mean(multispectral_matrix, 3);
centered_multispectral = multispectral_matrix - repmat(mean_multispectral, 1, 1, size(multispectral_matrix, 3));

% Reshape the centered_multispectral to a 2D matrix for PCA
reshaped_multispectral = reshape(centered_multispectral, [], size(centered_multispectral, 3));

% Compute the covariance matrix using X^T*X
covariance_matrix_multispectral = (reshaped_multispectral' * reshaped_multispectral) / (size(reshaped_multispectral, 1));

% Perform PCA on the covariance matrix
[eigenvectors_multispectral, eigenvalues_multispectral] = eig(covariance_matrix_multispectral);

% Sort eigenvectors based on eigenvalues in descending order
[eigenvalues_multispectral, sorted_indices_multispectral] = sort(diag(eigenvalues_multispectral), 'descend');
eigenvectors_multispectral = eigenvectors_multispectral(:, sorted_indices_multispectral);

% Number of principal components to keep
num_components = 4;

% Use the first 'num_components' principal components for pansharpening
top_pcs_multispectral = eigenvectors_multispectral(:, 1:num_components);

% Replace the top component of the multispectral image to be put through
% inverse PCA
reshaped_multispectral(:, 1) = reshape(panchromatic_image, m*n, 1);

% Pansharpening using the top principal components
pansharpened_components = reshape(reshaped_multispectral * top_pcs_multispectral, m, n, num_components) + mean_multispectral;

% Rescale the original images
NIR = rescale(multispectral_images{1, 1}(:,:,4));
R = rescale(multispectral_images{1, 1}(:,:,3));
G = rescale(multispectral_images{1, 1}(:,:,2));
B = rescale(multispectral_images{1, 1}(:,:,1));

% Calculate the NDVI for the original image
NDVI = ((NIR) - (R)) ./ ((NIR) + (R));

% Rescale the pansharpened images
NIR_P = rescale(pansharpened_components(:,:,4));
R_P = rescale(pansharpened_components(:,:,3));
G_P = rescale(pansharpened_components(:,:,2));
B_P = rescale(pansharpened_components(:,:,1));

% Calculate the NDVI for the pansharpened image
NDVI_P =  (((NIR_P - R_P) ./ (NIR_P + R_P)));

% Rescale the images and convert them to grayscale
original_image = im2gray(uint8(rescale(cat(3, R, G, B), 0, 255)));
pca_pansharpened_image = im2gray(uint8(rescale(cat(3, R_P, G_P, B_P), 0, 255)));

% -------------------------------------------------------------------------

% Calculate the denominator for the brovey method
multispectral_bands = multispectral_images{1,1};
denominator = multispectral_bands(:,:,1) + multispectral_bands(:,:,2) + multispectral_bands(:,:,3) + multispectral_bands(:,:,4);

% Calculate the bands for the brovey method
fused_band1 = (multispectral_bands(:,:,1).*panchromatic_image)./denominator;
fused_band2 = (multispectral_bands(:,:,2).*panchromatic_image)./denominator;
fused_band3 = (multispectral_bands(:,:,3).*panchromatic_image)./denominator;
fused_band4 = (multispectral_bands(:,:,4).*panchromatic_image)./denominator;

% Rescale the image and convert it to grayscale
fused_image = im2gray(uint8(rescale(cat(3, fused_band1, fused_band2, fused_band3), 0, 255)));

% Calculate the NDVI for the brovey image
brovey_NDVI = ((fused_band4) - (fused_band3)) ./ ((fused_band4) + (fused_band3));

% -------------------------------------------------------------------------

% Create pansharpened, brovey, and original images
pansharpened_image = cat(3, R_P, G_P, B_P);
brovey_image = cat(3, fused_band1, fused_band2, fused_band3);
starting_image = cat(3, R, G, B);

% Create pansharpened, brovey, and original images with all four bands
pansharpened_image_full = cat(3, R_P, G_P, B_P, NIR_P);
brovey_image_full = cat(3, fused_band1, fused_band2, fused_band3, fused_band4);
starting_image_full = cat(3, R, G, B, NIR);

% Calculate the Euclidian distance for the PCA image
for i=1:4  
    pca_euclidian_distance(i) = sqrt(sum((imhist(starting_image_full(:,:,i)) - imhist(pansharpened_image_full(:,:,i))).^2));
end

% Calculate the Euclidian distance for the Brovey image
for i=1:4  
    bro_euclidian_distance(i) = sqrt(sum((imhist(starting_image_full(:,:,i)) - imhist(brovey_image_full(:,:,i))).^2));
end

PCA_SAM_Matrix = [];
Brovey_SAM_Matrix = [];

% Calculate the spectral angle matrix for the PCA Image
for i=1:length(starting_image(:,:,:))
    for j=1:length(starting_image(:,:,:))
        starting_slice = reshape(starting_image(i,j,:), 3, 1);
        pansharpened_slice = reshape(pansharpened_image(i,j,:), 3, 1);
        pan_cos_insert = (pansharpened_slice' * starting_slice) / (norm(pansharpened_slice) * norm(starting_slice));
        PCA_SAM_Matrix(i, j) = acosd(pan_cos_insert);
    end
end

% Calculate the spectral angle matrix for the Brovey Image
for i=1:length(starting_image(:,:,:))
    for j=1:length(starting_image(:,:,:))
        starting_slice = reshape(starting_image(i,j,:), 3, 1);
        brovey_slice = reshape(brovey_image(i,j,:), 3, 1);
        bro_cos_insert = (brovey_slice' * starting_slice) / (norm(brovey_slice) * norm(starting_slice));
        Brovey_SAM_Matrix(i, j) = acosd(bro_cos_insert);
    end
end

% Reshape the euclidian distance matrices
pca_euclidian_distance_reshaped = [pca_euclidian_distance(1), pca_euclidian_distance(2), pca_euclidian_distance(3), pca_euclidian_distance(4)]';
bro_euclidian_distance_reshaped = [bro_euclidian_distance(1), bro_euclidian_distance(2), bro_euclidian_distance(3), bro_euclidian_distance(4)]';

% Calculate the average of the spectral angle matrices
PCA_SAM_Matrix_mean = mean(mean(PCA_SAM_Matrix));
Brovey_SAM_Matrix_mean = mean(mean(Brovey_SAM_Matrix));

% Create the images labels
images = {"Red"; "Green"; "Blue"; "NIR"};

% Create the tables for the spectral angle averages
spatial_table = table(images, pca_euclidian_distance_reshaped, bro_euclidian_distance_reshaped);
spectral_table = table(PCA_SAM_Matrix_mean, Brovey_SAM_Matrix_mean);

% -------------------------------------------------------------------------

figure;
subplot(2, 4, 1);
imshow(cat(3, R, zeros(m, n), zeros(m, n)));
title(['Original Red']);
subplot(2, 4, 2);
imshow(cat(3, zeros(m, n), G, zeros(m, n)));
title(['Original Green']);
subplot(2, 4, 3);
imshow(cat(3, zeros(m, n), zeros(m, n), B));
title(['Original Blue']);
subplot(2, 4, 4);
imshow(cat(3, NIR, NIR, NIR));
title(['Original NIR']);
subplot(2, 4, 5);
imshow(cat(3, R_P, zeros(m, n), zeros(m, n)));
title(['PCA Red']);
subplot(2, 4, 6);
imshow(cat(3, zeros(m, n), G_P, zeros(m, n)));
title(['PCA Green']);
subplot(2, 4, 7);
imshow(cat(3, zeros(m, n), zeros(m, n), B_P));
title(['PCA Blue']);
subplot(2, 4, 8);
imshow(cat(3, NIR_P, NIR_P, NIR_P), []);
title(['PCA NIR']);

figure;
subplot(1, 4, 1);
imshow(original_image, []);
title(['Multispectral Image']);
subplot(1, 4, 2);
imshow(panchromatic_image, []);
title(['Panchromatic Image']);
subplot(1, 4, 3);
imshow(pca_pansharpened_image, []);
title(['PCA Pansharpened Image']);
subplot(1, 4, 4);
imshow(fused_image);
title(['Brovey Method Pansharpened Image']);

figure;
subplot(1, 3, 1);
imshow(NDVI);
title(['Mulitspectral NDVI'])
subplot(1, 3, 2);
imshow(NDVI_P, []);
title(['PCA Pansharpening NDVI']);
subplot(1, 3, 3);
imshow(brovey_NDVI);
title(['Brovey Pansharpening NDVI'])

figure;
subplot(1, 2, 1);
imshow(PCA_SAM_Matrix, []);
title(["PCA Spectral Similarity"]);
subtitle(num2str(PCA_SAM_Matrix_mean));
subplot(1, 2, 2);
imshow(Brovey_SAM_Matrix, []);
title(["Brovey Spectral Similarity"]);
subtitle(num2str(Brovey_SAM_Matrix_mean));

uitable(uifigure, 'Data', spatial_table, 'ColumnName', {'Image Slice', 'PCA Euclidian Distance', 'Brovey Euclidian Distance'}, "Position",[20 20 450 117]);
uitable(uifigure, 'Data', spectral_table, 'ColumnName', {'PCA Average Spectral Angle', 'Brovey Average Spectral Angle'}, "Position",[20 20 400 51]);
