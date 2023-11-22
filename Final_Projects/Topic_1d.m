
% Load the image dataset
load ms_IKONOS.mat;

% Load the image dataset
load pan_IKONOS.mat;

ms_IKONOS = struct2cell(load("ms_IKONOS.mat"));
pan_IKONOS = struct2cell(load("pan_IKONOS.mat"));

multispectral_images = cellfun(@double , ms_IKONOS, 'UniformOutput', false);
panchromatic_image = pan_IKONOS{1, 1};

% Resize multispectral images to match panchromatic image size
[m, n, channels, frames] = size(panchromatic_image);
multispectral_images = cellfun(@(x) imresize(x, [m, n]), multispectral_images, 'UniformOutput', false);

% Combine multispectral images into a single matrix
multispectral_matrix = cell2mat(reshape(multispectral_images, 1, 1, []));

%X is an L×N data matrix with L bands and N pixels;
mean_multispectral = mean(multispectral_matrix')';

% estimate the covariance matrix;
covariance_multispectral = multispectral_matrix*multispectral_matrix'/(m*n) - mean_multispectral*mean_multispectral'; 

% eigen-decomposition;
[e1 e2 e3] = eig(covariance_multispectral); 

% first PC
pc1 = e1(:, 4)'*(multispectral_matrix-mean_multispectral*ones(1,(m*n)));

% second PC
pc2 = e2(:, 4-1)'*(multispectral_matrix-mean_multispectral*ones(1,(m*n)));

% reshape for display
pc1 = reshape(pc1, [m, n]);

% display
figure; imshow(pc1,[]); 

