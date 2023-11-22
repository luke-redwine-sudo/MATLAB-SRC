close all

% Load the image dataset
load Sentinel2-data.mat;

Sentinel2_data = struct2cell(load("Sentinel2-data.mat"));

numSlices = 5;

chosenImage = 5;


% Load the image dataset
load ms_IKONOS.mat;

ms_IKONOS = load("ms_IKONOS.mat");

% Load the image dataset
load pan_IKONOS.mat;

pan_IKONOS = struct2cell(load("pan_IKONOS.mat"));

% Combine the multispectral and panchromatic data
ms_data = struct2cell(load("ms_IKONOS.mat"));

%ms_data = ms_data{1, 1}; % Assuming that the multispectral data is stored in a field named 'data'
pan_data = load("pan_IKONOS.mat"); % Assuming that the panchromatic data is stored in a field named 'data'

% Determine the layout based on the number of images
layout_rows = ceil(sqrt(numSlices));
layout_columns = ceil(numSlices / layout_rows);

% Display all images simultaneously
% figure;

% for slice = 1 : numSlices
%     subplot(layout_rows, layout_columns, slice);
%     imshow(Sentinel2_data{slice}, []);
%     title(['Image ' num2str(slice)]);
%     disp(['Displaying Image ' num2str(slice)]); % Add this line for debugging
% end

% Assume the first image is the panchromatic image, and the rest are multispectral
%panchromatic_image = double(Sentinel2_data{chosenImage});
%multispectral_images = cellfun(@double, Sentinel2_data([1:4, 6:end]), 'UniformOutput', false);
%multispectral_images = cellfun(@double, Sentinel2_data(2:end), 'UniformOutput', false);

multispectral_images = cellfun(@double , ms_data, 'UniformOutput', false);
panchromatic_image = pan_IKONOS{1, 1};

% Resize multispectral images to match panchromatic image size
[m, n] = size(panchromatic_image);
multispectral_images = cellfun(@(x) imresize(x, [m, n]), multispectral_images, 'UniformOutput', false);

% Combine multispectral images into a single matrix
multispectral_matrix = cell2mat(reshape(multispectral_images, 1, 1, []));

% Perform PCA on the multispectral images
mean_multispectral = mean(multispectral_matrix, 3);
centered_multispectral = multispectral_matrix - repmat(mean_multispectral, 1, 1, size(multispectral_matrix, 3));

% Reshape the centered_multispectral to a 2D matrix for PCA
reshaped_multispectral = reshape(centered_multispectral, [], size(centered_multispectral, 3));

reshaped_multispectral(:, 1) = reshape(panchromatic_image, m*n, 1);

% Compute the covariance matrix using X^T*X
covariance_matrix_multispectral = (reshaped_multispectral' * reshaped_multispectral) / (size(reshaped_multispectral, 1) - 1);

% Perform PCA on the covariance matrix
[eigenvectors_multispectral, eigenvalues_multispectral] = eig(covariance_matrix_multispectral);

% Sort eigenvectors based on eigenvalues in descending order
[eigenvalues_multispectral, sorted_indices_multispectral] = sort(diag(eigenvalues_multispectral), 'descend');
eigenvectors_multispectral = eigenvectors_multispectral(:, sorted_indices_multispectral);

% Number of principal components to keep
num_components = 3;

% Use the first 'num_components' principal components for pansharpening
top_pcs_multispectral = eigenvectors_multispectral(:, 1:num_components);

% Pansharpening using the top principal components
pansharpened_components = reshape(reshaped_multispectral * top_pcs_multispectral, m, n, num_components);
pansharpened_image = sum(pansharpened_components, 3) + panchromatic_image;

% Display the original and pansharpened images
figure;
subplot(1, 2, 1);
imshow(panchromatic_image, []);
title(['Panchromatic Image']);

subplot(1, 2, 2);
imshow(pansharpened_image, []);
title(['Pansharpened Image - Top ' num2str(num_components) ' PCs']);

R = rescale(multispectral_images{1, 1}(:,:,3));
G = rescale(multispectral_images{1, 1}(:,:,2));
B = rescale(multispectral_images{1, 1}(:,:,1));

figure;
subplot(1, 2, 1);
imshow(cat(3, R, G, B));
title(['Multispectral Color Image']);
subplot(1, 2, 2);
imshow(cat(3, pansharpened_components(:,:,1), pansharpened_components(:,:,2), pansharpened_components(:,:,3)));
title(['Pansharpened Image']);
