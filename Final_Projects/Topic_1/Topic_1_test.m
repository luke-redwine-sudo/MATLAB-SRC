paths = ['pca_pansharp:', 'common:'];
addpath(paths);

load ms_IKONOS.mat
A = ms(:,:, 1:3);
A = double(A)/255;
%A = cellfun(@(x) imresize(x, [m, n]), ms_data, 'UniformOutput', false);


load pan_IKONOS.mat
B = double(pan)/255;

sharped = solve_pansharp(A, B);

figure, imshow(imresize(A, 2));
figure, imshow(sharped);

rmpath(paths);