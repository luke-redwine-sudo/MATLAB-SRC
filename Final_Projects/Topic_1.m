% Load the image dataset
load Sentinel2-data.mat;

image_data_structure = struct2cell(load("Sentinel2-data.mat"));

% Initialize dimensions
rows = 521;
columns = 508;
numSlices = 13;

image_list = cell(1, numSlices);

% Extract image slices from the loaded structure and resize to a common size
common_size = [rows, columns];
for slice = 1 : numSlices
    image = double(image_data_structure{slice, 1});
    resized_image = imresize(image, common_size);
    image_list{slice} = resized_image;
end

% Reshape the resized images into column vectors
image_matrix = reshape(cell2mat(image_list), [], numSlices);

% Subtract the mean from each dimension
mean_image = mean(image_matrix, 2);
image_matrix_centered = image_matrix - mean_image;

% Compute the covariance matrix using X^T*X
covariance_matrix = (image_matrix_centered' * image_matrix_centered) / (numSlices - 1);

% Compute the eigenvectors and eigenvalues
[eigenvectors, eigenvalues] = eig(covariance_matrix);

% Sort eigenvectors based on eigenvalues in descending order
[eigenvalues, sorted_indices] = sort(diag(eigenvalues), 'descend');
eigenvectors = eigenvectors(:, sorted_indices);

% Number of principal components to keep (you can adjust this)
num_components = 3;

% Use the first 'num_components' principal components
top_pcs = eigenvectors(:, 1:num_components);

% Transform the top principal components back into images
reconstructed_images = (image_matrix_centered * top_pcs) * top_pcs' + mean_image;

% Display the reconstructed images in separate figures
for slice = 1 : num_components
    figure;
    subplot(1, 2, 1);
    imshow(reshape(image_matrix(:, slice), [rows, columns]), []);
    title(['Original Image ' num2str(slice)]);
    
    subplot(1, 2, 2);
    imshow(reshape(reconstructed_images(:, slice), [rows, columns]), []);
    title(['Reconstructed Image ' num2str(slice) ' - Top ' num2str(num_components) ' PCs']);
end

% Determine the layout based on the number of images
layout_rows = ceil(sqrt(numSlices));
layout_columns = ceil(numSlices / layout_rows);

% Display all images simultaneously
figure;

for slice = 1 : numSlices
    subplot(layout_rows, layout_columns, slice);
    imshow(image_list{slice}, []);
    title(['Image ' num2str(slice)]);
    disp(['Displaying Image ' num2str(slice)]); % Add this line for debugging
end