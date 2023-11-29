
addpath('C:\Users\redwi\Documents\Final_Projects');
ncdisp("GEO2.nc")

imagesc(flip(rot90(uint8(h5read("C:\Users\redwi\Documents\Final_Projects\GEO2.nc", "/QA")), 3), 2));