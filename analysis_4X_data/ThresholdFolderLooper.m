%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loops through a folder structure, where at the bottom each folder contains
% 8-bit images to be thresholded using the histogram decomposition functions. This script
% assumes that all images are at an equal depth in the folder structure and is only written
% for convenience. The folder structure we assume is:
% Master directory:
%     Directories with treatments (e.g. drugs)
%         Directories with replicates
%             Images to be processed
% Your folder structure might be different and each level might correspond
% to a different thing.
%
% This script outputs a structure called output but also stores a text file with its
% values in the master directory.
%
% Author: Jose L. Cadavid, University of Toronto, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
% Initialize output
output = struct([]);
% Master index
ii = 1;

% Get current (master) folder: This folder has subfolders (drugs), each
% with subfolders (reps), each with N number of image files to be processed
master = uigetdir([], 'Please select the master folder to be processed');
disp(strcat('Processing folder: ', master));

% Create text file (must add a meaningful name, and this function should
% loop across folders
fid = fopen(strcat(master,'\Results_MATLAB.txt'),'w');
% Print headers - add date
fprintf(fid, "%s \n", strcat("Results of image segmentation performed on ",master," on ",date));

fprintf(fid,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s \n", 'Treatment', 'Rep', 'Image',...
                'TreatmentID', 'RepID', 'ImgID',...
                'Percentage_coverage', 'Optimal_threshold', 'Flag', 'Percentage_saturated');

% List of files in master folder
listM = dir(fullfile(master,'*'));
% List of subfolders (drugs) in master folder
drugs = setdiff({listM([listM.isdir]).name},{'.','..'});

% Loop through drugs
for drg = 1:numel(drugs)
    disp(strcat('Processing treatment: ', drugs{drg}));
    % Get list of files in drug folder
    listD = dir(fullfile(master,drugs{drg},'*'));
    % List of subfolders in drug (reps)
    reps = setdiff({listD([listD.isdir]).name},{'.','..'});

    % Loop through reps
    for rp = 1:numel(reps)
        disp(strcat('Processing rep: ', reps{rp}));
        % Get list of files in rep
        listR = dir(fullfile(master,drugs{drg},reps{rp},'*'));
        % Get list of actual images in rep (images are NOT folders. Folders must ONLY contain images)
        imgs = {listR(~[listR.isdir]).name};

        % Loop through images
        for im = 1:numel(imgs)
            output(ii).drug = drugs{drg};
            output(ii).rep = reps{rp};
            output(ii).imgName = erase(imgs{im},'.tif');

            % Include drug, rep and image counters for easy averaging and
            % sorting
            output(ii).drugN = drg;
            output(ii).repN = rp;
            output(ii).imgN = im;

            % Threshold image

            % Read image
            img = imread(fullfile(master,drugs{drg},reps{rp},imgs{im}),'tif');
            % Operate
            [pArea, tOpt,flag,~] = thresholdImage(img,0);
            % Store
            output(ii).flag = flag;
            output(ii).pArea = pArea;
            output(ii).tOpt = tOpt;
            % Percentage of saturated pixels
            output(ii).pSat = sum(img==255,'all')/numel(img)*100;

            % Write in text file
            fprintf(fid,"%s\t %s\t %s\t %u\t %u\t %u\t %.2f\t %0.2f\t %u\t %.2f \n", ...
                drugs{drg}, reps{rp}, imgs{im},...
                drg, rp, im, pArea, tOpt, flag, output(ii).pSat);
            % Update master index
            ii = ii+1;
        end
    end
end

% Close text file
fclose(fid);
disp('Done!');
