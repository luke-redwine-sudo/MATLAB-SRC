clc; clear all; close all
%% Loading File
pic0 = imread('Mona_Lisa-256x256.jpeg');     % Read image
picgrey = rgb2gray(pic0);     % Coverting to grayscale
pic = double(picgrey);        % Converting from uint8 to double for processing
%% Compression Ratios
DivBy = [16];  % Section Square Size/Compression Ratio
M = length(pic(:,1)); N = length(pic(1,:)); %   Row by Col
SNR = zeros(length(DivBy), 2);
%% Compression Loop
for k = 1:length(DivBy)
% Reshaping File
N_fix = fix(N/DivBy(k))*DivBy(k) % How many rows/cols to keep
pic = double(picgrey(1:N_fix,1:N_fix));
% Initializing    
np = DivBy(k)*DivBy(k);    % Number of pixels in each section
delta = DivBy(k) - 1;      % The Jump
L = (N_fix*N_fix)/(DivBy(k)^2); % Number of sections
l = 0;
% Matrix initializing
matrix = zeros(DivBy(k),DivBy(k),L); 
CompPic = zeros(DivBy(k),DivBy(k),L);
CompPic2 = zeros(N_fix,N_fix);
%% PCA Compression
for i = 1:DivBy(k):(N_fix-delta)
    for j = 1:DivBy(k):(N_fix-delta)
    l = l+1;
    matrix(:,:,l) = pic(i:i+delta,j:j+delta); 
    end
end
z = reshape(matrix, [np l]);
Avg = mean(z,2);
Z = z-Avg;
Covar = (Z*Z')/N ;
[V, D] = eig(Covar);
y = V(:,(length(V(1,:))-DivBy(k)+1):end)'*Z;

%% PCA Reconstruction
new = V(:,(length(V(1,:))-DivBy(k)+1):end)*y+Avg;
j = 0;
DivBy(k)
for j = 1:length(new(1,:))
   CompPic(:,:,j) = reshape(new(:,j), [DivBy(k) DivBy(k)]);
end
l = 0;
for i = 1:DivBy(k):(N_fix-delta)
   for j = 1:DivBy(k):(N_fix-delta)
    l = l+1;
    CompPic2(i:i+delta,j:j+delta) = abs(CompPic(:,:,l));
    end
end
%% PCA SNR
%% 
fhat = sum(sum(CompPic2));
f = sum(sum(pic));
SNR(k,1) = (fhat^2)/(fhat-f)^2; 

figure(k)
imshow(uint8(CompPic2),[])
title(sprintf('Reconstructed %d:1',DivBy(k)))

%% DCT
% Compress
N_fix
DCTPic = dct2(picgrey,N_fix,N_fix);
% Removing Negatives
DCTPic(abs(DCTPic) < 10) = 0;
% Reconstruct
fix(N/DivBy(k))
DCTPic2 = idct2(DCTPic,fix(N/DivBy(k)),fix(N/DivBy(k)));
DCTPic2 = rescale(DCTPic2);
figure(k+4)
imshow(DCTPic2, [])
title(sprintf('Reconstructed DCT of %1.0f:1', DivBy(k)))
%% DCT SNR
fhat = sum(sum(255*double(DCTPic2)));
SNR(k,2) = ((fhat)^2)/((fhat-f)^2)




end

