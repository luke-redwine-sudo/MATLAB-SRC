clear
close all

x = imread('bird.jpg');
x = reshape(x,[131*175 3]);
x = double(x);

x = x';

m = mean(x')';

v = x*x'/22925 - m*m';

%correlation coefficient
coff = zeros(3);

for i=1:3
    for j=1:3
        coff(i,j) = v(i,j)/sqrt(v(i,j))/sqrt(v(j,j));
    end
end

% eigenvector martix
[v1 v2] = eig(v);

v1*v1'
v1'*v1

% mutual orthogonal
v1(:,1)'*v1(:,2)

pca1 = v1(:,3)'*(x-m*ones(1,22925));
pca2 = v1(:,2)'*(x-m*ones(1,22925));
pca3 = v1(:,1)'*(x-m*ones(1,22925));

pca1 = reshape(pca1, [131 175]);
pca2 = reshape(pca2, [131 175]);
pca3 = reshape(pca3, [131 175]);

subplot(1,3,1); imshow(pca1, []);
subplot(1,3,2); imshow(pca2, []);
subplot(1,3,3); imshow(pca3, []);

color_image = zeros(131, 175);
color_image(:,:,1) = pca1;
color_image(:,:,2) = pca2;
color_image(:,:,3) = pca3;

figure
imshow(color_image, [])
