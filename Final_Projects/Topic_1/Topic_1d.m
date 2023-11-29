
% Load the image dataset
load ms_IKONOS.mat;

% Load the image dataset
load pan_IKONOS.mat;

ms_IKONOS = struct2cell(load("ms_IKONOS.mat"));
pan_IKONOS = struct2cell(load("pan_IKONOS.mat"));

ms_IKONOS2 = struct2cell(load("imagery\cropped\1_ms.mat"));
pan_IKONOS = struct2cell(load("imagery\cropped\1_pan.mat"));

multispectral_images = ms_IKONOS{1, 1};
panchromatic_image = pan_IKONOS{1, 1};

% Resize multispectral images to match panchromatic image size
[m, n, l] = size(panchromatic_image);

multispectral_images = reshape(multispectral_images, m*n, []);

%X is an L×N data matrix with L bands and N pixels;
mean_multispectral = mean(multispectral_images')';

% estimate the covariance matrix;
cov1 = multispectral_images*multispectral_images'/(m*n);
covariance_multispectral = cov1 - mean_multispectral*mean_multispectral'; 

% eigen-decomposition;
[e1 e2 e3] = eig(covariance_multispectral); 
%% 

% first PC
pc1a = multispectral_images-mean_multispectral*ones(1,(m*n));
pc1 = e1(:, l)'*pc1a;

% second PC
%pc2 = e2(:, 4-1)'*(multispectral_images-mean_multispectral*ones(1,(m*n)));

% reshape for display
pc1 = reshape(pc1, [m, n]);

% display
figure; imshow(pc1,[]); 

