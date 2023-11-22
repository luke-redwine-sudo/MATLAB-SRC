a = imread('dark.jpg'); % read an image
figure; image(a); % display a color image
a = double(a);
a = (a(:,:,1)+a(:,:,2)+a(:,:,3))/3;

[N,X] = hist(a(:), 255);
figure; stem(X, N)

figure; imhist(a/max(a(:)))
[Count, XX] = imhist(a/max(a(:)));
figure; stem(XX, Count)

a1 = round(a);
p = zeros(1,256);
for i = 1:131
for j = 1:175
    p(a1(i,j)) = p(a1(i,j)) + 1;
end
end
figure; stem(1:256, p)