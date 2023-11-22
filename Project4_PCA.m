%% ECE 8473 Digital Image Processing: Project 4, PCA
% By Rebecca Garcia
% Dated: October 24th, 2022
clc; clear all;
%% Loading File
pic0 = imread('woman_hat.jpg');     % Read image
pic = double(pic0);        % Converting from uint8 to double for processing
%% Reshaping File 
DivBy = 4;  % Section Square Size/Compression Ratio
M = length(pic(:,1)); N = length(pic(1,:)); %   Row by Col
M_fix = fix(M/DivBy)*DivBy; % How many rows to keep
N_fix = fix(N/DivBy)*DivBy; % How many cols to keep
pic = pic(1:M_fix,1:N_fix);
%% PCA
n = 0; m = 0;       % Initializing Indexes
CompPic = zeros(fix(M/DivBy),fix(N/DivBy));
N = DivBy*DivBy;    % Number of pixels in each section
delta = DivBy - 1;  % The Jump
% Loop
for i = 1:DivBy:(M_fix-delta)
    m = m+1;
    for j = 1:DivBy:(N_fix-delta)
        n = n+1;
        z = reshape(pic(i:i+delta,j:j+delta), [N 1]);   % Reshaping 
        Avg = mean(z);      % Finding Average
        Z = z-Avg;          % Removing Average
        Covar = (Z*Z')/N ;  % Covariance Matrix
        [V, D] = eig(Covar);
        [D] = diag(D);                  % EigenValues 
        v = V(:,end)         % First Eigen Vector
        CompPic(m,n) = abs(v'*Z);                   
    end 
    n = 0;
end
figure
hold on
imshow(uint8(CompPic),[])


