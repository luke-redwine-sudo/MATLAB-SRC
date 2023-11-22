% Load the image dataset
load moffett_AVIRIS.mat;

image_data_structure = struct2cell(load("moffett_AVIRIS.mat"));

% Initialize dimensions
rows = 521;
columns = 508;
numSlices = 13;

image_data = struct('data', zeros(rows, columns, numSlices));

% Extract image slices from the loaded structure
for slice = 1 : 2
    image = uint8(image_data_structure{slice, 1});
    imshow(image(1:rows, 1:columns));
end



% Extract image slices from the loaded structure
for slice = 1 : numSlices
    slice_name = sprintf('B%d', slice); % Assuming slices are named 'slice1', 'slice2', ...
    var = image_data_structure{slice, 1}(rows, columns, numSlices)
end
