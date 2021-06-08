%Adapted from: https://www.mathworks.com/matlabcentral/answers/437494-how-to-loop-through-all-files-in-subfolders-in-a-main-folder
clear all
clc
%Initialize output
output=struct([]);
%Master index
ii=1;

%Get current (master) folder: This folder has subfolders (drugs), each
%with subfolders (reps), each with N number of image files to be processed



master=cd;
%List of files in master folder
list_m=dir(fullfile(master,'*'));
%List of subfolders (drugs) in master folder
drugs=setdiff({list_m([list_m.isdir]).name},{'.','..'});

%Loop through drugs
for i=1:numel(drugs)
    %Get list of files in drug
    list_d=dir(fullfile(master,drugs{i},'*'));
    %List of subfolders in drug (reps)
    reps=setdiff({list_d([list_d.isdir]).name},{'.','..'});
    
    %Loop through reps
    for j=1:numel(reps)
        %Get list of files in rep
        list_r=dir(fullfile(master,drugs{i},reps{j},'*'));
        %Get list of actual images in rep (images are NOT folders. Folders must ONLY contain images)
        imgs={list_r(~[list_r.isdir]).name};
        
        %Loop through images
        for k=1:numel(imgs)
            output(ii).drug=drugs{i};
            output(ii).rep=reps{j};
            output(ii).im_name=erase(imgs{k},'.tif');
            
            %Include drug, rep and image counters for easy averaging and
            %sorting
            output(ii).drugN=i;
            output(ii).repN=j;
            output(ii).imgN=k;
            
            %Threshold image
            
            %Read image
            im=imread(fullfile(master,drugs{i},reps{j},imgs{k}),'tif');
            %Operate
            [pArea, t_opt,flag,~]=thresholdImage(im,0);
            %Store
            output(ii).flag=flag;
            output(ii).pArea=pArea;
            output(ii).tOpt=t_opt;
            
            output(ii).pSat=sum(im==255,'all')/numel(im)*100; %percentage of saturated pixels
            %Update master index
            ii=ii+1;
        end
        
    end
    
end



%%%%% After data is sorted use this to average across experiments (average
%%%%% technical replicates

%  i=1;
% avg_data=[];
% for i=1:11
% subM=DATA(DATA(:,1)==i,:); %Submatrix of drug
% %Average for each exp
% Nexp=max(subM(:,2));
% for j=1:Nexp
% avg_exp=mean(subM(subM(:,2)==j,4:6),1);
% N_imgs=max(subM(subM(:,2)==j,3));
% %Store
% avg_data(ii,:)=[i,j,N_imgs,avg_exp];
% %update
% ii=ii+1;
% end
% end