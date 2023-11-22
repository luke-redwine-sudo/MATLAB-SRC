a = imread('book.tif'); figure; imshow(a,[]);
a = double(a);
ya = fft2(a);
figure; imshow(abs(ya),[])
figure; imshow(log(abs(ya)),[])

a2 = a;
for i = 1:size(a,1)
    for j = 1:size(a,2)
        a2(i,j) = a(i,j)*(-1)^(i+j);
    end
end
ya2 = fft2(a2);
figure; imshow(abs(ya2),[])
figure; imshow(log(abs(ya2)),[])

ya3=fftshift(fft2(a)); figure; imshow(log(abs(ya3)),[])
