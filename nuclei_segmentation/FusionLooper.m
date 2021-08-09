%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loops through a folder that contains resampled 3D
% single-channel images (SAA, GFP, Nuclei) of different samples (identified
% by a number on the image).
%
% This script outputs a structure called output but also stores a text file with its
% values in the master directory selected.
%
% Author: Jose L. Cadavid, University of Toronto, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% Master index
ii = 1;

% Get folder with images to be analyzed
master = uigetdir([], 'Please select the master folder to be processed');
master = strcat(master, filesep);
disp(strcat('Processing folder: ', master));
% Get text file exported by FIJI with image characteristics
[filename, pathname] = uigetfile({'*.txt'},'Select text file with info');
% Read table
imgTable = readtable(strcat(pathname,filename), 'Delimiter', {'\t'});
% Number of images is the height of the table
nImgs = height(imgTable);

% Loop through images

for i = 1:nImgs
    disp(strcat('Processing image #',num2str(i)));
    % Load channels
    imgSAA = readTiffStack(strcat(master,'SAA_reslice',num2str(i),'.tif'));
    imgNuc = readTiffStack(strcat(master,'Nuc_reslice',num2str(i),'.tif'));
    imgGFP = readTiffStack(strcat(master,'GFP_reslice',num2str(i),'.tif'));
    % Process
    scale = imgTable.Height_um_(i);
    [~, nNuclei, nGFP, nSAA, nGFPSAA] = getNuclearFractions(imgNuc, imgGFP, imgSAA, scale);
     % Store in table
    imgTable.nNuclei(i) = nNuclei;
    imgTable.nGFP(i) = nGFP;
    imgTable.nSAA(i) = nSAA;
    imgTable.nGFPSAA(i) = nGFPSAA;
end

% Save table
writetable(imgTable, strcat(pathname,erase(filename, '.txt'),'_processed.txt'), 'Delimiter', '\t');
disp('Done!');
